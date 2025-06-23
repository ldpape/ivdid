cap program drop sivdid
program define sivdid, eclass 
syntax  [if] [in] [,   Y(string) D(string) Z(string) first(real 0)  median  periods(real 1) arg_cohort_treatment_date(string) arg_time(string) exponential controls(varlist) keep(real 5)  ]        
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
tempvar total_dates  sum_dates nb_dates sum_nb_dates
qui: bys `arg_cohort_treatment_date' `arg_time'  : gen `total_dates' = (_n==`keep') if `_sample' == 1 // & (`arg_cohort_treatment_date' > 0) // at least ten observations per cohort and date 
//qui: replace `_sample' = 0 if `total_dates' == 0 & `_sample' == 1 & (`arg_cohort_treatment_date' > 0) // at least ten observations per cohort and date 
 qui: bys `arg_cohort_treatment_date'  : egen `sum_dates' = sum(`total_dates')   if `_sample' == 1 //& (`arg_cohort_treatment_date' > 0) // check if condition validates across all dates 
 qui: bys `arg_cohort_treatment_date' `arg_time'  : gen `nb_dates' = (_n==1) if `_sample' == 1 //& (`arg_cohort_treatment_date' > 0) // count number of dates
 qui: bys `arg_cohort_treatment_date'  : egen `sum_nb_dates' = sum(`nb_dates')   if `_sample' == 1 //& (`arg_cohort_treatment_date' > 0) // check if condition validates across all dates 
 qui : replace `_sample' = 0  if (`sum_dates'< `sum_nb_dates')  & `_sample'==1 //& (`arg_cohort_treatment_date' > 0) // keep only valid cohorts 
/*   BEGIN ESTIMATION      */
tempvar beta_hat
cap: gen `beta_hat' = .
local cow = 0
qui: levelsof `arg_cohort_treatment_date' if `arg_cohort_treatment_date'>0 & `_sample'==1 // loop through cohorts and dates
local cohort_groups `r(levels)'
local nb_cohorts `r(r)'
foreach t_group in `cohort_groups' {
local cow = `cow' + 1
		forvalues t_ell = 0(1)`periods' {
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
cap: scalar Ferror_ivreg2 = (e(F) < `first') // | (abs(e(b)[1,1])<1e-8)  | (abs(e(b)[1,2])<1e-8)  | (abs(e(b)[1,3])<1e-8) | (abs(e(b)[1,4])<1e-8) 
scalar Ferror_immediate = _rc 		
cap: scalar Ferror_ivreg = (abs(e(b)[1,1]) == .) | (abs(e(b)[1,2])==.)  | (abs(e(b)[1,3])==.) | (abs(e(b)[1,4])==.)  
    * select estimation method (linear or ppml)					
if "`exponential'" == "" {
		* LINEAR MODEL (OLS/IV)
cap:    ivreg2 `y' ( `d' = `z') `Variable_G'  `Variable_T' `controls'  if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample'  
scalar error_immediate = _rc 
						}
					else {
		* IV POISSON OPTION
cap:  ivreg2 `lnY' ( `d' = `z') `Variable_G'  `Variable_T' `controls'  if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample'
// cap: noisily:  ivpoisson  gmm `y' `Variable_G' `Variable_T' `controls' (`d'= `z') if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample' , from(e(b))  multiplicative  technique(nr)  onestep  conv_maxiter(20) 
 cap:  ivpois `y' `Variable_G'  `Variable_T' `controls' if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample', endog(`d') exog(`z') from(e(b))
 scalar error_immediate = _rc 
						 }		
     * extract treatment effect and potential errors in second stage 
//  cap: scalar beta_iv_`cow'_`t_ell' = _b[`d']  
 cap: replace `beta_hat' = _b[`d'] if ( (`arg_cohort_treatment_date'  == (`t_group')) ) & (`arg_time' == (`t_group' + `t_ell')) & `_sample'
//  cap: scalar error_ivreg2 = (abs(e(b)[1,1]) < 1e-40) | (abs(e(b)[1,2])<1e-40)  | (abs(e(b)[1,3])<1e-40) | (abs(e(b)[1,4])<1e-40) 
 cap: scalar error_ivreg = (abs(e(b)[1,1]) == .) | (abs(e(b)[1,2])==.)  | (abs(e(b)[1,3])==.) | (abs(e(b)[1,4])==.)  
	 * remove corrupted estimates and associated sample 
 if error_ivreg2 == 1 | error_ivreg == 1 | error_immediate>0 |  Ferror_ivreg2 == 1 | Ferror_ivreg == 1 | Ferror_immediate>0 {  
//  	di "Could not estimate the folowing group-time treament effect:"
//  	di "Group-Time: " `cow' "-" `t_ell'
    cap: replace `beta_hat' = . if ( (`arg_cohort_treatment_date'  == (`t_group')) ) & (`arg_time' == (`t_group' + `t_ell')) & `_sample'
	qui: replace `_sample' = 0 if (`arg_cohort_treatment_date' == `t_group') & (`arg_time' == (`t_group' + `t_ell')) & `_sample'==1
			}
 				}
					}	 			
/*   REPORT AVERAGE TREATMENT EFFECTS      */
** report to eclass 
ereturn clear 
qui: sum `_sample' if `_sample'
scalar obs_actual = r(N)
qui: sum `copy' if `copy'
scalar N = r(N)
cap : drop sample 
cap : gen sample = `_sample'
ereturn post, esample(`copy') obs(`=N') depname(`y')
** report in console 
di  "********************************************************************"
di  "*                 SUMMARY OF ESTIMATED EFFECTS                     *"
di  "********************************************************************"

if "`exponential'" == "exponential" { 
di "-- Estimates based on endogenous count model (ivpois) -- "
}
di "Dynamic Average Effects :"
scalar beta_error = 0
* report per period effect
forvalues t_ell = 0(1)`periods' {
qui: sum `beta_hat' if `_sample' & (`arg_cohort_treatment_date' != 0) & (`arg_time' == (`arg_cohort_treatment_date' + `t_ell' )) , de 
 if (r(N) > 0) {
	if "`median'" == "median"{
	di "Median Treatment Effect for Period `t_ell' : " r(p50)
	ereturn scalar beta_median_`t_ell' = r(p50)
			}
	else {
	di "Mean Treatment Effect for Period : " r(mean)
	}
	ereturn scalar beta_mean_`t_ell' = r(mean)	
	}
	else {
		di "BETA_`t_ell' is not identified"
		if "`median'" == "median"{
	ereturn scalar beta_median_`t_ell' = .
			}
		else{
	di "BETA_`t_ell' is not identified"	
		}
	ereturn scalar beta_mean`t_ell' = .	
	}
}
* report overall effect 
qui: sum `beta_hat' if `_sample' & (`arg_cohort_treatment_date' != 0) , de
 if (r(N) > 0) {
	if "`median'" == "median"{
	di "Overall Median Treatment Effect : " r(p50)
	ereturn scalar beta_overall_median = r(p50)
	di "Overall Mean Treatment Effect : " r(mean // provide even if median option
			}
	else {
	di "Overall Mean Treatment Effect : " r(mean)
	}
	ereturn scalar beta_overall_mean = r(mean)		
			  }
	else {	
		if "`median'" == "median"{
		di "Beta_median is not identified"
	ereturn scalar beta_overall_median = .
			}
		else{
		di "BETA_mean is not identified"
		}
	ereturn scalar beta_overall_mean = .
	scalar beta_error = 1
	}
* ereturn (eclass) estimates as scalars (to do: matrix form)
ereturn scalar error = beta_error
ereturn scalar N_obs = obs_actual
ereturn scalar keep = `keep'
ereturn scalar first = `first'
ereturn scalar periods = `periods'
di "Number of Observations: " obs_actual
end
