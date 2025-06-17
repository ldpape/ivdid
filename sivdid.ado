cap program drop miyaji
program define miyaji, eclass 
syntax varlist [if] [in] [aweight pweight fweight iweight] [,   Y(string) D(string) Z(string) first(real 1)   periods(real 1) arg_cohort_treatment_date(string) arg_time(string) expo(string) controls(varlist) Keep(real 5)  ]        
/*         PARSE TEXT       */
marksample _sample
	set buildfvinfo on
	tempvar copy 
	qui : gen `copy' = `_sample' // used for jackknife to identify initial sample
/*   CHECK SAMPLE VALIDITY  */
tempvar lnY
qui: gen `lnY' = asinh(`y') // used for exponential option 
	* keep sample with dates corresponding to request
qui: replace `_sample' = 0 if ((`arg_time' > (`arg_cohort_treatment_date' + `periods')) | (`arg_time' < (`arg_cohort_treatment_date' - 1))) & (`arg_cohort_treatment_date' > 0)
tempvar total_dates  sum_dates 
qui: bys `arg_cohort_treatment_date' `arg_time'  : gen `total_dates' = (_n==`keep') if `_sample' == 1 & (`arg_cohort_treatment_date' > 0) // at least ten observations per cohort and date 
qui: bys `arg_cohort_treatment_date'  : egen `sum_dates' = sum(`total_dates')   if `_sample' == 1 & (`arg_cohort_treatment_date' > 0) // check if condition validates across all dates 
qui : replace `_sample' = 0  if `sum_dates'< (1+`periods')  & `_sample'==1 & (`arg_cohort_treatment_date' > 0) // keep only valid cohorts 
/*   BEGIN ESTIMATION      */
local cow = 0
qui: levelsof `arg_cohort_treatment_date' if `arg_cohort_treatment_date'>0 & `_sample'==1 // loop through cohorts and dates
local cohort_groups `r(levels)'
local nb_cohorts `r(r)'
foreach t_group in `cohort_groups' {
local cow = `cow' + 1
		forvalues t_ell = 0(1)`periods' {
	* create an empty treatment effect scalar 
scalar beta_iv_`cow'_`t_ell' = .
	* scalars to record potential identification failure 
scalar error_ivreg2 = 0
scalar error_ivreg = 0
scalar error_immediate = 0
scalar Ferror_ivreg2 = 0
scalar Ferror_ivreg = 0
scalar Ferror_immediate = 0
	* create dummies for time and cohort 
qui: tempvar Variable_G 
qui: tempvar Variable_T
qui: gen `Variable_G' = (`arg_cohort_treatment_date'  == (`t_group'))	if `_sample'
qui: gen `Variable_T' = (`arg_time' == (`t_group' + `t_ell')) if `_sample'
    * check first-stage : record error if first-stage (t< X) is too weak or not enough variation for fixed effects
cap:  reg  `d' `z' `Variable_G'  `Variable_T' `controls'  if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample'  
cap: scalar Ferror_ivreg2 = (e(F) < `first') | (abs(e(b)[1,2])<1e-8)  | (abs(e(b)[1,3])<1e-8) | (abs(e(b)[1,4])<1e-8) 
scalar Ferror_immediate = _rc 		
cap: scalar Ferror_ivreg = (abs(e(b)[1,1]) == .) | (abs(e(b)[1,2])==.)  | (abs(e(b)[1,3])==.) | (abs(e(b)[1,4])==.)  
    * select estimation method (linear or ppml)					
if "`expo'" != "exponential" {
		* LINEAR MODEL (OLS/IV)
cap:    ivreg2 `y' ( `d' = `z') `Variable_G'  `Variable_T' `controls'  if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample'  
scalar error_immediate = _rc 		
						}
					else {
		* IV POISSON OPTION
cap:  ivreg2 `lnY' ( `d' = `z') `Variable_G'  `Variable_T' `controls'  if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample'
cap:  ivpoisson  gmm `y' `Variable_G' `Variable_T' `controls' (`d'= `z') if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample' , from(e(b))  multiplicative  technique(nr)  onestep  conv_maxiter(300) 
 scalar error_immediate = _rc 
						 }		
     * extract treatment effect and potential errors in second stage 
 cap: scalar beta_iv_`cow'_`t_ell' = _b[`d']
 cap: scalar error_ivreg2 = (abs(e(b)[1,1]) < 1e-8) | (abs(e(b)[1,2])<1e-8)  | (abs(e(b)[1,3])<1e-8) | (abs(e(b)[1,4])<1e-8) 
 cap: scalar error_ivreg = (abs(e(b)[1,1]) == .) | (abs(e(b)[1,2])==.)  | (abs(e(b)[1,3])==.) | (abs(e(b)[1,4])==.)  
	 * remove corrupted estimates and associated sample 
 if error_ivreg2 == 1 | error_ivreg == 1 | error_immediate>0 |  Ferror_ivreg2 == 1 | Ferror_ivreg == 1 | Ferror_immediate>0 {  
//  	di "Could not estimate the folowing group-time treament effect:"
//  	di "Group-Time: " `cow' "-" `t_ell'
	scalar beta_iv_`cow'_`t_ell' = 0
	qui: replace `_sample' = 0 if (`arg_cohort_treatment_date' == `t_group') & (`arg_time' == (`t_group' + `t_ell')) & `_sample'==1
			}
 				}
					}	 
/*   CALCULATE WEIGHTS BY TREATMENT EFFECTS      */
local cow = 0
foreach t_group in `cohort_groups' {
local cow = `cow' + 1
		forvalues t_ell = 0(1)`periods' {
* calculate number of observation treated for date and cohort 
qui: sum `d' if (`arg_cohort_treatment_date' == (`t_group')) & (`arg_time' == (`t_group' + `t_ell')  )  & `_sample' == 1
qui: scalar share = r(N)
* calculate number of observation for date
qui: sum `d' if (`arg_cohort_treatment_date' != 0) & (`arg_time' == (`arg_cohort_treatment_date' + `t_ell' ))   & `_sample' == 1
qui: scalar total = r(N) 
* calculate weight by time : weight_iv_cohort_time
if total == 0 { // no identified 
qui: scalar weight_iv_`cow'_`t_ell' = 0
}
else {
qui: scalar weight_iv_`cow'_`t_ell' = share/total  
}
* calculate weight for overall sample  : total_weight_iv_cohort_time
qui: sum `arg_cohort_treatment_date' if `arg_cohort_treatment_date' != 0 & (`arg_time' >= `arg_cohort_treatment_date') & (`arg_time'<= (`arg_cohort_treatment_date' + `periods'))  & `_sample' == 1
qui: scalar total_sample = r(N) 
if total == 0 {
qui: scalar total_weight_iv_`cow'_`t_ell' =  0
}
else {
qui: scalar total_weight_iv_`cow'_`t_ell' =  share/total_sample 
}
 				}
					}		
/*   CALCULATE AVERAGE TREATMENT EFFECTS      */
scalar beta_AVG = 0 // preparation scalars
scalar total_weight = 0
scalar count_ = 0
scalar sum_ = 0
* iterate across time and groups 
		forvalues t_ell = 0(1)`periods' {
scalar	beta_`t_ell' = 0
scalar  weight_`t_ell' = 0
local   cow = 0
				foreach t_group in `cohort_groups' {
local cow = `cow' + 1
scalar beta_`t_ell' = beta_`t_ell' + beta_iv_`cow'_`t_ell'*weight_iv_`cow'_`t_ell' // add weights 
scalar beta_AVG = beta_AVG +  beta_iv_`cow'_`t_ell'*total_weight_iv_`cow'_`t_ell' // overall treatment effect
scalar weight_`t_ell' =  weight_`t_ell' + weight_iv_`cow'_`t_ell' // used to check correct implemenation 
scalar total_weight =  total_weight + total_weight_iv_`cow'_`t_ell' // used to check correct implemenation 
 ** calculate weights by variance of period  
  if beta_iv_`cow'_`t_ell' != 0{ 
  	 scalar count_ = count_ + 1 
 	 scalar sum_ = sum_ + beta_iv_`cow'_`t_ell'
  }
		}
			}
/*   REPORT AVERAGE TREATMENT EFFECTS      */
** report to eclass 
ereturn clear 
qui: sum `_sample' if `_sample'
scalar N = r(N)
ereturn post, esample(`_sample') obs(`=N') depname(`y')
** report in console 
di in red "********************************************************************"
di in red "*                 SUMMARY OF ESTIMATED EFFECTS                     *"
di in red "********************************************************************"

if "`exponential'" == "exponential" { 
di in red "-- Estimates based on endogenous count model (ivpois) -- "
}
di in red "Dynamic Average Effects :"
	scalar beta_error = 0
forvalues t_ell = 0(1)`periods' {
if beta_`t_ell' == 0 { // return error if no beta_ell
di "BETA_`t_ell' is not identified" 
	scalar beta_`t_ell' = . // for jackknife/bootstrap, report identification failure
//	scalar beta_error = 1
	}
	else {		
	di "BETA_`t_ell' = " beta_`t_ell'
	}
ereturn scalar beta_`t_ell' = beta_`t_ell'
ereturn scalar weight_`t_ell' = weight_`t_ell'
}
if beta_AVG == 0 { // return error or parameters 
	di "BETA_AVG is not identified"
	scalar beta_AVG = .
	scalar beta_error = 1
}
else {
di in red "Average Effect per Period : "
	di "BETA_AVG =" beta_AVG	
di in red "Unweighted Average Effect per Period : "
	di "BETA_UWAVG =" sum_/count_	
}
* ereturn (eclass) estimates as scalars (to do: matrix form)
ereturn scalar beta_AVG = beta_AVG
ereturn scalar total_weight = total_weight
ereturn scalar beta_UWAVG = sum_/count_
ereturn scalar error = beta_error
end
