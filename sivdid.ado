cap program drop sivdid
program define sivdid, eclass 
syntax [if] [in] [,   Y(string) D(string) Z(string) first(real 0)  event_study  periods(real 0) cohort_treatment_date(string) time(string) exponential controls(varlist) keep(real 0) permanent graph_options(string) ]        
/*         PARSE TEXT       */
marksample _sample 
tempvar copy 
// qui : gen `copy' = `_sample' // used for jackknife to identify initial sample
* drop missing obs 
foreach x of varlist `y' `d' `z' `cohort_treatment_date' `time' {
qui: replace `_sample' = 0 if missing(`x')
}	
/*   CHECK SAMPLE VALIDITY  */
	* for treated, drop periods outside of cohort + periods, or periods - 1
qui: replace `_sample' = 0 if ((`time' > (`cohort_treatment_date' + `periods')) | (`time' < (`cohort_treatment_date' - 1))) & (`cohort_treatment_date' > 0) & `_sample'
	* remove observations when there are no treated cohorts (useful for control group)
 tempvar  max_obs
 qui: bys `time' : egen `max_obs' = sum((`_sample'==1)*(`cohort_treatment_date'>0)) 
 qui replace `_sample' = 0 if `max_obs' == 0 & `_sample'
	* keep consistent set of cohorts across all _sample 
if "`permanent'" != ""{
	* check if there are at least `keep' amount of observations per cohort at each relevant moments in time
	* then, check if this condition is valided across moments in time
	* drop a cohort which is not available for the whole length of time
tempvar total_dates  sum_dates nb_dates sum_nb_dates
 qui: bys  `_sample' `cohort_treatment_date' `time'  : gen `total_dates' = (_n==`keep') if `_sample' == 1
 qui: bys  `_sample' `cohort_treatment_date'         : egen `sum_dates' = sum(`total_dates')   if `_sample' == 1 
 qui: bys  `_sample' `cohort_treatment_date' `time'  : gen `nb_dates' = (_n==1) if `_sample' == 1 
 qui: bys  `_sample' `cohort_treatment_date'         : egen `sum_nb_dates' = sum(`nb_dates')   if `_sample' == 1 
 qui : replace `_sample' = 0  if (`sum_dates'< `sum_nb_dates')  & `_sample'==1 
}
	* gen variable for exponential 
	if "`exponential'"!=""{
tempvar lnY
qui: gen `lnY' = asinh(`y') if `_sample' 
	}
/*   BEGIN ESTIMATION      */
tempvar beta_hat se_hat
cap: gen `se_hat' = .
cap: gen `beta_hat' = .
local cow = 0
qui: levelsof `cohort_treatment_date' if `cohort_treatment_date'>0 & `_sample' // loop through cohorts and dates
local cohort_groups `r(levels)'
if "`cohort_groups'" == ""{
	di in red "Error: no valid cohorts identified post data cleaning."
}
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
scalar FS = 0
	* create dummies for time and cohort 
qui: tempvar Variable_G 
qui: tempvar Variable_T
qui: gen `Variable_G' = (`cohort_treatment_date'  == (`t_group'))	if `_sample'
qui: gen `Variable_T' = (`time' == (`t_group' + `t_ell')) if `_sample'
    * check first-stage : record error if first-stage (t< X) is too weak or not enough variation for fixed effects
cap: reg  `d' `z' `Variable_G'  `Variable_T' `controls'  if ( (`cohort_treatment_date'  == (`t_group')) | (`cohort_treatment_date' == 0)) & (`time' == (`t_group' + `t_ell')   | (`time'== `t_group' - 1)) & `_sample'  
cap: scalar Ferror_ivreg2 = (e(F) < `first') // | (abs(e(b)[1,1])<1e-6)  | (abs(e(b)[1,2])<1e-6)  | (abs(e(b)[1,3])<1e-6) | (abs(e(b)[1,4])<1e-6) 
scalar Ferror_immediate = _rc 		
cap: scalar FS = _b[`z']
cap: scalar Ferror_ivreg = (abs(e(b)[1,1]) == .) | (abs(e(b)[1,2])==.)  | (abs(e(b)[1,3])==.) | (abs(e(b)[1,4])==.)  
    * select estimation method (linear or ppml)					
if "`exponential'" == "" {
		* LINEAR MODEL (OLS/IV)
cap: ivreg2 `y' ( `d' = `z') `Variable_G'  `Variable_T' `controls'  if ( (`cohort_treatment_date'  == (`t_group')) | (`cohort_treatment_date' == 0)) & (`time' == (`t_group' + `t_ell')   | (`time'== `t_group' - 1)) & `_sample'  
scalar error_immediate = _rc 
cap: scalar SS = _b[`d']
						}
					else {
		* IV POISSON OPTION
cap:  ivreg2 `lnY' ( `d' = `z') `Variable_G'  `Variable_T' `controls'  if ( (`cohort_treatment_date'  == (`t_group')) | (`cohort_treatment_date' == 0)) & (`time' == (`t_group' + `t_ell')   | (`time'== `t_group' - 1)) & `_sample'
// cap: noisily:  ivpoisson  gmm `y' `Variable_G' `Variable_T' `controls' (`d'= `z') if ( (`cohort_treatment_date'  == (`t_group')) | (`cohort_treatment_date' == 0)) & (`time' == (`t_group' + `t_ell')   | (`time'== `t_group' - 1)) & `_sample' , from(e(b))  multiplicative  technique(nr)  onestep  conv_maxiter(20) 
 cap:  ivpois `y' `Variable_G'  `Variable_T' `controls' if ( (`cohort_treatment_date'  == (`t_group')) | (`cohort_treatment_date' == 0)) & (`time' == (`t_group' + `t_ell')   | (`time'== `t_group' - 1)) & `_sample', endog(`d') exog(`z') from(e(b))
 scalar error_immediate = _rc 
						 }		
     * extract treatment effect and potential errors in second stage 
cap: replace `beta_hat' = _b[`d'] if ( (`cohort_treatment_date'  == (`t_group')) ) & (`time' == (`t_group' + `t_ell')) & `_sample'
// cap: scalar error_ivreg2 = (abs(e(b)[1,1]) < 1e-6) | (abs(e(b)[1,2])<1e-6)  | (abs(e(b)[1,3])<1e-6) | (abs(e(b)[1,4])<1e-6) 
cap: scalar error_ivreg = (abs(e(b)[1,1]) == .) | (abs(e(b)[1,2])==.)  | (abs(e(b)[1,3])==.) | (abs(e(b)[1,4])==.)  
	 * remove corrupted estimates and associated sample 
 if error_ivreg2 == 1 | error_ivreg == 1 | error_immediate>0 |  Ferror_ivreg2 == 1 | Ferror_ivreg == 1 | Ferror_immediate>0 | FS==0 {  
//  	di "Could not estimate the folowing group-time treament effect:"
//  	di "Group-Time: " `cow' "-" `t_ell'
    cap: replace `beta_hat' = . if ( (`cohort_treatment_date'  == (`t_group')) ) & (`time' == (`t_group' + `t_ell')) & `_sample'
	qui: replace `_sample' = 0 if (`cohort_treatment_date' == `t_group') & (`time' == (`t_group' + `t_ell')) & `_sample'
	scalar SS = . 
	scalar FS = .
			}
 				
/*   CALCULATE STANDARD ERRORS OF EACH CLATT     */
if "`exponential'" == ""  { // option only valid for linear 2SLS 
	*tempvars 
tempvar  delta_`cow'_`t_ell'  phi_`cow'_`t_ell' phi_sqr_`cow'_`t_ell'  se_`cow'_`t_ell'
	*calculate share of treated
qui: sum `_sample' if `_sample'  
scalar total_N = r(N)
qui: sum  `Variable_G' if  (`time' == (`t_group' + `t_ell')) & ( (`cohort_treatment_date'  == (`t_group')) | (`cohort_treatment_date' == 0)) & `_sample'  
scalar E_T1 = r(mean)*r(N)/total_N // share treated post 
scalar E_C1 = (1-r(mean))*r(N)/total_N  // share control post 
qui: sum  `Variable_G' if  (`time' == (`t_group' - 1)) & ( (`cohort_treatment_date'  == (`t_group')) | (`cohort_treatment_date' == 0)) & `_sample'  
scalar E_T0 = r(mean)*r(N)/total_N // share treated pre 
scalar E_C0 = (1-r(mean))*r(N)/total_N  // share control pre
	* calculate residuals
qui: gen `delta_`cow'_`t_ell'' = `y' -  SS*`d' if  `_sample'
	* calculate conditional means 
qui:sum `delta_`cow'_`t_ell''  if (`cohort_treatment_date'  == `t_group') &(`time' == (`t_group' + `t_ell')) & `_sample' 
scalar D_T1 = r(mean) // resid for treated post 
qui:sum `delta_`cow'_`t_ell''  if (`cohort_treatment_date' == 0) &(`time' == (`t_group' + `t_ell')) & `_sample' 
scalar D_C1 = r(mean) // resid for control post
qui: sum `delta_`cow'_`t_ell''  if (`cohort_treatment_date'  == `t_group') &(`time' == (`t_group' -1)) & `_sample' 
scalar D_T0 = r(mean) // resid for treated pre 
qui: sum `delta_`cow'_`t_ell''  if (`cohort_treatment_date' == 0) &(`time' == (`t_group' -1)) & `_sample' 
scalar D_C0 = r(mean) // resid for control pre
	* calculate influence function phi 
qui: gen `phi_`cow'_`t_ell'' = (((`cohort_treatment_date'  == (`t_group')) | (`cohort_treatment_date' == 0)) & (`time' == (`t_group' + `t_ell')   | (`time'==  (`t_group' - 1))) ) *	(1/FS)*(   (`cohort_treatment_date'  == (`t_group'))*(`time' == (`t_group' + `t_ell'))*( `delta_`cow'_`t_ell''-D_T1)/E_T1 + (`cohort_treatment_date'  == (`t_group'))*(`time' == (`t_group' -1))*( `delta_`cow'_`t_ell''-D_T0)/E_T0 + (`cohort_treatment_date'  == 0)*(`time' == (`t_group' + `t_ell'))*( `delta_`cow'_`t_ell''-D_C1)/E_C1 - (`cohort_treatment_date'  == (0))*(`time' == (`t_group' -1))*( `delta_`cow'_`t_ell''-D_C0)/E_C0)  if `_sample'
	* calculate standard error of each clatt
qui: gen `phi_sqr_`cow'_`t_ell'' = (`phi_`cow'_`t_ell'')^2  if `_sample'
qui: sum `phi_sqr_`cow'_`t_ell'' if  `_sample'
qui: gen `se_`cow'_`t_ell'' = sqrt(r(mean)/r(N))
cap: replace `se_hat' = sqrt(r(mean)/r(N)) if ( (`cohort_treatment_date'  == (`t_group')) ) & (`time' == (`t_group' + `t_ell')) & `_sample'
//  sum se_`cow'_`t_ell'
// sum `beta_hat' if `_sample' &  (`time' == (`t_group'+ `t_ell')) & ( (`cohort_treatment_date'  == (`t_group')))
}									
	}
		}
if "`exponential'" == ""   { // option only valid for linear 2SLS 
/*   CALCULATE STANDARD ERRORS OF EACH TIME PERIOD    */
forvalues t_ell = 0(1)`periods' {
tempvar se_theta_`t_ell' influence_`t_ell'  influence_sqr_`t_ell' 
qui: sum `beta_hat' if `_sample' & (`cohort_treatment_date' != 0) & (`time' == (`cohort_treatment_date' + `t_ell' ) )
scalar total_`t_ell' = r(N)
qui: gen `influence_`t_ell'' = 0 if `_sample'
local cow = 0
foreach t_group in `cohort_groups' {
local cow = `cow'  + 1
tempvar psi_`cow'_`t_ell' 
qui: gen `psi_`cow'_`t_ell'' = ((`cohort_treatment_date'  == (`t_group')) & (`time' == (`cohort_treatment_date' + `t_ell' ) ))  if `_sample' 
qui: sum `psi_`cow'_`t_ell'' if `_sample' & (`cohort_treatment_date'!= 0 ) & (`time' == (`cohort_treatment_date' + `t_ell' ) )
qui: replace `psi_`cow'_`t_ell'' = `psi_`cow'_`t_ell'' - r(mean) if `_sample' & (`cohort_treatment_date' != 0) & (`time' == (`cohort_treatment_date' + `t_ell' ) )
qui: sum `beta_hat' if  (`cohort_treatment_date'  == (`t_group')) & (`time' == (`t_group' + `t_ell')) & `_sample'  
if r(N)>0{
qui: replace `influence_`t_ell'' = `influence_`t_ell''  + ((r(N)/total_`t_ell')*(`phi_`cow'_`t_ell'') + r(mean)*`psi_`cow'_`t_ell'')  if `_sample'  
         }
	}
qui: gen `influence_sqr_`t_ell''= (`influence_`t_ell'')^2  if `_sample'
qui: sum `influence_sqr_`t_ell'' if `_sample'
qui: gen `se_theta_`t_ell''= sqrt(r(mean)/r(N))
//sum se_theta_`t_ell'
}

/*   CALCULATE STANDARD ERRORS FOR FINAL THETA  */
qui: sum `beta_hat' if `_sample' & (`cohort_treatment_date' != 0) & (`time' >= (`cohort_treatment_date')) 
scalar total_final = r(N)
tempvar influence_final  influence_sqr_final se_theta_final
qui: gen  `influence_final' = 0
forvalues t_ell = 0(1)`periods' {
local cow = 0
foreach t_group in `cohort_groups' {
local cow = `cow'  + 1
qui: sum `beta_hat' if  (`cohort_treatment_date'  == (`t_group')) & (`time' == (`t_group' + `t_ell')) & `_sample'  
if r(N)>0{
qui: replace `influence_final' = `influence_final'  + ((r(N)/total_final)*(`phi_`cow'_`t_ell'') + r(mean)*`psi_`cow'_`t_ell'')  if  `_sample'
		 }	
	}
}
qui: gen `influence_sqr_final' = (`influence_final')^2 if `_sample'
qui: sum `influence_sqr_final' if `_sample'
qui: gen `se_theta_final' = sqrt(r(mean)/r(N))
												}
/*   REPORT AVERAGE TREATMENT EFFECTS      */
** report to eclass 
ereturn clear 
qui: sum `_sample' if `_sample'
scalar obs_actual = r(N)
qui: gen `copy' = `_sample'
// scalar N = r(N)
//  cap : drop sample 
//  cap : gen sample = `_sample' // actual sample used, made available directly
ereturn post, esample(`copy') obs(`=obs_actual') depname(`y')
** report in console 
di  "----------------------------------------------------"
di  "                   LATE SUMMARY                     "
di  "----------------------------------------------------"
scalar beta_error = 0
scalar Ncoef = `periods'+2
matrix b = J(Ncoef,1,.)
matrix se = J(Ncoef,1,.)
* Fill period effects 
forvalues t_ell = 0(1)`periods' {
qui: sum `beta_hat' if `_sample' & (`cohort_treatment_date' != 0) & (`time' == (`cohort_treatment_date' + `t_ell' )) 
scalar beta_error = beta_error + r(N)==0
matrix b[`t_ell'+1, 1] = r(mean)
qui: sum  `se_theta_`t_ell'' 
matrix se[`t_ell'+1, 1] = r(mean)
}
* Fill with avg. effect (theta)
qui: sum `beta_hat' if `_sample' & (`cohort_treatment_date' != 0) & (`time' >= (`cohort_treatment_date')) 
matrix b[Ncoef, 1] = r(mean)
qui: sum `se_theta_final'
matrix se[Ncoef, 1] = r(mean)

if "`exponential'" != ""  { // not calculated for exponential model
matrix se = J(Ncoef,1,.)
}
* Prepare to export 
tempvar print_beta_hat print_se_hat lower_ci upper_ci ell
cap: gen `print_beta_hat' = .
cap: gen `print_se_hat' = .
cap: gen `ell' = .
local names 
di as text "----------------------------------------------------"
di as text "     LATE     |   Coef.    Std. Err.     z     P>|z|"
di as text "----------------------------------------------------"
forvalues t_ell = 0(1)`periods' {
        local coef = b[`t_ell'+1,1]
        local ster = se[`t_ell'+1,1]
        local ratio = `coef'/`ster'
        local p = 2 * (1 - normal(abs(`ratio')))
        local varname = "Period `t_ell'" 
	    local names `names' LATE_`t_ell'
        di as res %10s "`varname'" ///
           "   " %9.3f `coef' ///
           "   " %9.3f `ster' ///
           "   " %6.2f `ratio' ///
           "   " %6.3f `p'
cap: replace `print_beta_hat' = b[`t_ell'+1,1] if _n == (`t_ell'+1)
cap: replace `print_se_hat' = se[`t_ell'+1,1] if _n == (`t_ell'+1)
cap: replace `ell' = (`t_ell') if _n == (`t_ell'+1)
    }
        local coef = b[Ncoef,1]
        local ster = se[Ncoef,1]
        local ratio = `coef'/`ster'
        local p = 2 * (1 - normal(abs(`ratio')))
        local varname = "Total Avg."
	    local names `names' LATE_AVG
        di as res %10s "`varname'" ///
           "   " %9.3f `coef' ///
           "   " %9.3f `ster' ///
           "   " %6.2f `ratio' ///
           "   " %6.3f `p'
    di as text "----------------------------------------------------"
matrix rownames b = `names'
matrix rownames se = `names'
ereturn matrix late = b
ereturn matrix std = se 
ereturn scalar error = beta_error
ereturn scalar keep = `keep'
ereturn scalar first = `first'
ereturn scalar periods = `periods'
cap : drop _late_hat
cap : gen _late_hat = `beta_hat'
cap : drop _late_hat_se
cap : gen _late_hat_se = `se_hat'
di "Number of Observations: " obs_actual
if "`event_study'" != ""{
cap : gen `lower_ci' = `print_beta_hat' - 1.96*`print_se_hat'
cap : gen `upper_ci' = `print_beta_hat' + 1.96*`print_se_hat'
tw (scatter `print_beta_hat' `ell', mcolor(eltblue) msize(2.5) yline(`coef', lpattern(dash) lwidth(0.35) lcolor(blue)) xlabel(,grid glcolor(gray) glpattern(dash) glwidth(.1))  ylabel( , grid glcolor(gray) glpattern(dash) glwidth(.1))  ) (rcap `upper_ci' `lower_ci' `ell' , lcolor(eltblue) lwidth(0.25) lpattern(solid)) , graphregion(color(white)) plotregion(color(white)) bgcolor(white)  xtitle("Periods Since Treatment")  ytitle("Estimated per Period Avg. LATE")  legend(off) yline(0, lpattern(solid) lwidth(0.35) lcolor(red)) `graph_options' 
						}
end
