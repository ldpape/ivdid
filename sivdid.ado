cap program drop sivdid
program define sivdid, eclass 
syntax [if] [in] [,   Y(string) D(string) Z(string) first(real 0)  median  periods(real 1) arg_cohort_treatment_date(string) arg_time(string) exponential controls(varlist) keep(real 1)  ]        
/*         PARSE TEXT       */
marksample _sample
// 	set buildfvinfo on
tempvar copy 
qui : gen `copy' = `_sample' // used for jackknife to identify initial sample
* drop missing obs 
foreach x of varlist `y' `d' `z' `arg_cohort_treatment_date' `arg_time' {
qui: replace `_sample' = 0 if missing(`x')
}	
/*   CHECK SAMPLE VALIDITY  */
	* keep sample with dates corresponding to request
qui: replace `_sample' = 0 if ((`arg_time' > (`arg_cohort_treatment_date' + `periods')) | (`arg_time' < (`arg_cohort_treatment_date' - 1))) & (`arg_cohort_treatment_date' > 0) & `_sample'
tempvar total_dates  sum_dates nb_dates sum_nb_dates
qui: bys `_sample' `arg_cohort_treatment_date' `arg_time'  : gen `total_dates' = (_n==`keep') if `_sample' // & (`arg_cohort_treatment_date' > 0) // at least ten observations per cohort and date 
 qui: bys `_sample' `arg_cohort_treatment_date'  : egen `sum_dates' = sum(`total_dates')   if `_sample' //& (`arg_cohort_treatment_date' > 0) // check if condition validates across all dates 
 qui: bys `_sample' `arg_cohort_treatment_date' `arg_time'  : gen `nb_dates' = (_n==1) if `_sample' //& (`arg_cohort_treatment_date' > 0) // count number of dates
 qui: bys `_sample' `arg_cohort_treatment_date'  : egen `sum_nb_dates' = sum(`nb_dates')   if `_sample' //& (`arg_cohort_treatment_date' > 0) // check if condition validates across all dates 
 qui : replace `_sample' = 0  if (`sum_dates'< `sum_nb_dates')  & `_sample' //& (`arg_cohort_treatment_date' > 0) // keep only valid cohorts 
 tempvar maxdate 
 qui: egen `maxdate' = max(`arg_cohort_treatment_date')
 qui : replace `_sample' = 0 if `arg_cohort_treatment_date'==0 &  ((`arg_time' > (`maxdate' + `periods'))) & `_sample'
tempvar lnY
qui: gen `lnY' = asinh(`y') if `_sample' // used for exponential option 
/*   BEGIN ESTIMATION      */
tempvar beta_hat
cap: gen `beta_hat' = .
local cow = 0
qui: levelsof `arg_cohort_treatment_date' if `arg_cohort_treatment_date'>0 & `_sample' // loop through cohorts and dates
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
qui: gen `Variable_G' = (`arg_cohort_treatment_date'  == (`t_group'))	if `_sample'
qui: gen `Variable_T' = (`arg_time' == (`t_group' + `t_ell')) if `_sample'
    * check first-stage : record error if first-stage (t< X) is too weak or not enough variation for fixed effects
cap: reg  `d' `z' `Variable_G'  `Variable_T' `controls'  if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample'  
cap: scalar Ferror_ivreg2 = (e(F) < `first') // | (abs(e(b)[1,1])<1e-6)  | (abs(e(b)[1,2])<1e-6)  | (abs(e(b)[1,3])<1e-6) | (abs(e(b)[1,4])<1e-6) 
scalar Ferror_immediate = _rc 		
cap: scalar FS = _b[`z']
cap: scalar Ferror_ivreg = (abs(e(b)[1,1]) == .) | (abs(e(b)[1,2])==.)  | (abs(e(b)[1,3])==.) | (abs(e(b)[1,4])==.)  
    * select estimation method (linear or ppml)					
if "`exponential'" == "" {
		* LINEAR MODEL (OLS/IV)
cap: ivreg2 `y' ( `d' = `z') `Variable_G'  `Variable_T' `controls'  if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample'  
scalar error_immediate = _rc 
cap: scalar SS = _b[`d']
						}
					else {
		* IV POISSON OPTION
cap:  ivreg2 `lnY' ( `d' = `z') `Variable_G'  `Variable_T' `controls'  if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample'
// cap: noisily:  ivpoisson  gmm `y' `Variable_G' `Variable_T' `controls' (`d'= `z') if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample' , from(e(b))  multiplicative  technique(nr)  onestep  conv_maxiter(20) 
 cap:  ivpois `y' `Variable_G'  `Variable_T' `controls' if ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'== `t_group' - 1)) & `_sample', endog(`d') exog(`z') from(e(b))
 scalar error_immediate = _rc 
						 }		
     * extract treatment effect and potential errors in second stage 
cap: replace `beta_hat' = _b[`d'] if ( (`arg_cohort_treatment_date'  == (`t_group')) ) & (`arg_time' == (`t_group' + `t_ell')) & `_sample'
// cap: scalar error_ivreg2 = (abs(e(b)[1,1]) < 1e-6) | (abs(e(b)[1,2])<1e-6)  | (abs(e(b)[1,3])<1e-6) | (abs(e(b)[1,4])<1e-6) 
cap: scalar error_ivreg = (abs(e(b)[1,1]) == .) | (abs(e(b)[1,2])==.)  | (abs(e(b)[1,3])==.) | (abs(e(b)[1,4])==.)  
	 * remove corrupted estimates and associated sample 
 if error_ivreg2 == 1 | error_ivreg == 1 | error_immediate>0 |  Ferror_ivreg2 == 1 | Ferror_ivreg == 1 | Ferror_immediate>0 | FS==0 {  
//  	di "Could not estimate the folowing group-time treament effect:"
//  	di "Group-Time: " `cow' "-" `t_ell'
    cap: replace `beta_hat' = . if ( (`arg_cohort_treatment_date'  == (`t_group')) ) & (`arg_time' == (`t_group' + `t_ell')) & `_sample'
	qui: replace `_sample' = 0 if (`arg_cohort_treatment_date' == `t_group') & (`arg_time' == (`t_group' + `t_ell')) & `_sample'
	scalar SS = . 
	scalar FS = .
			}
 				
/*   CALCULATE STANDARD ERRORS OF EACH CLATT     */
if "`exponential'" == "" & "`median'" == ""  { // option only valid for linear 2SLS 
	*tempvars 
tempvar  delta_`cow'_`t_ell'  phi_`cow'_`t_ell' phi_sqr_`cow'_`t_ell'  se_`cow'_`t_ell'
	*calculate share of treated
qui: sum `_sample' if `_sample'  
scalar total_N = r(N)
qui: sum  `Variable_G' if  (`arg_time' == (`t_group' + `t_ell')) & ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & `_sample'  
scalar E_T1 = r(mean)*r(N)/total_N // share treated post 
scalar E_C1 = (1-r(mean))*r(N)/total_N  // share control post 
qui: sum  `Variable_G' if  (`arg_time' == (`t_group' - 1)) & ( (`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & `_sample'  
scalar E_T0 = r(mean)*r(N)/total_N // share treated pre 
scalar E_C0 = (1-r(mean))*r(N)/total_N  // share control pre
	* calculate residuals
qui: gen delta_`cow'_`t_ell' = `y' -  SS*`d' if  `_sample'
	* calculate conditional means 
qui:sum delta_`cow'_`t_ell'  if (`arg_cohort_treatment_date'  == `t_group') &(`arg_time' == (`t_group' + `t_ell')) & `_sample' 
scalar D_T1 = r(mean) // resid for treated post 
qui:sum delta_`cow'_`t_ell'  if (`arg_cohort_treatment_date' == 0) &(`arg_time' == (`t_group' + `t_ell')) & `_sample' 
scalar D_C1 = r(mean) // resid for control post
qui: sum delta_`cow'_`t_ell'  if (`arg_cohort_treatment_date'  == `t_group') &(`arg_time' == (`t_group' -1)) & `_sample' 
scalar D_T0 = r(mean) // resid for treated pre 
qui: sum delta_`cow'_`t_ell'  if (`arg_cohort_treatment_date' == 0) &(`arg_time' == (`t_group' -1)) & `_sample' 
scalar D_C0 = r(mean) // resid for control pre
	* calculate influence function phi 
qui: gen phi_`cow'_`t_ell' = (((`arg_cohort_treatment_date'  == (`t_group')) | (`arg_cohort_treatment_date' == 0)) & (`arg_time' == (`t_group' + `t_ell')   | (`arg_time'==  (`t_group' - 1))) ) *	(1/FS)*(   (`arg_cohort_treatment_date'  == (`t_group'))*(`arg_time' == (`t_group' + `t_ell'))*( delta_`cow'_`t_ell'-D_T1)/E_T1 + (`arg_cohort_treatment_date'  == (`t_group'))*(`arg_time' == (`t_group' -1))*( delta_`cow'_`t_ell'-D_T0)/E_T0 + (`arg_cohort_treatment_date'  == 0)*(`arg_time' == (`t_group' + `t_ell'))*( delta_`cow'_`t_ell'-D_C1)/E_C1 - (`arg_cohort_treatment_date'  == (0))*(`arg_time' == (`t_group' -1))*( delta_`cow'_`t_ell'-D_C0)/E_C0)  if `_sample'
	* calculate standard error of each clatt
qui: gen phi_sqr_`cow'_`t_ell' = (phi_`cow'_`t_ell')^2  if `_sample'
qui: sum phi_sqr_`cow'_`t_ell' if  `_sample'
qui: gen se_`cow'_`t_ell' = sqrt(r(mean)/r(N))
// sum se_`cow'_`t_ell'
// sum `beta_hat' if `_sample' &  (`arg_time' == (`t_group'+ `t_ell')) & ( (`arg_cohort_treatment_date'  == (`t_group')))
}									
	}
		}
/*   CALCULATE STANDARD ERRORS OF EACH TIME PERIOD    */
forvalues t_ell = 0(1)`periods' {
tempvar se_theta_`t_ell'
qui: sum `beta_hat' if `_sample' & (`arg_cohort_treatment_date' != 0) & (`arg_time' == (`arg_cohort_treatment_date' + `t_ell' ) )
scalar total_`t_ell' = r(N)
tempvar influence_`t_ell'  influence_sqr_`t_ell' 
qui: gen  influence_`t_ell' = 0 if `_sample'
local cow = 0
foreach t_group in `cohort_groups' {
local cow = `cow'  + 1
tempvar psi_`cow'_`t_ell' 
qui: gen psi_`cow'_`t_ell' = ((`arg_cohort_treatment_date'  == (`t_group')) & (`arg_time' == (`arg_cohort_treatment_date' + `t_ell' ) ))  if `_sample' 
qui: sum psi_`cow'_`t_ell' if `_sample' & (`arg_cohort_treatment_date'!= 0 ) & (`arg_time' == (`arg_cohort_treatment_date' + `t_ell' ) )
qui: replace psi_`cow'_`t_ell' = psi_`cow'_`t_ell' - r(mean) if `_sample' & (`arg_cohort_treatment_date' != 0) & (`arg_time' == (`arg_cohort_treatment_date' + `t_ell' ) )
qui: sum `beta_hat' if  (`arg_cohort_treatment_date'  == (`t_group')) & (`arg_time' == (`t_group' + `t_ell')) & `_sample'  
if r(N)>0{
qui: replace influence_`t_ell' = influence_`t_ell'  + ((r(N)/total_`t_ell')*(phi_`cow'_`t_ell') + r(mean)*psi_`cow'_`t_ell')  if `_sample'  
         }
	}
qui: gen influence_sqr_`t_ell' = (influence_`t_ell')^2  if `_sample'
qui: sum influence_sqr_`t_ell' if `_sample'
qui: gen se_theta_`t_ell' = sqrt(r(mean)/r(N))
//sum se_theta_`t_ell'
}

/*   CALCULATE STANDARD ERRORS FOR FINAL THETA  */
qui: sum `beta_hat' if `_sample' & (`arg_cohort_treatment_date' != 0) & (`arg_time' >= (`arg_cohort_treatment_date')) 
scalar total_final = r(N)
tempvar influence_final  influence_sqr_final se_theta_final
qui: gen  influence_final = 0
forvalues t_ell = 0(1)`periods' {
local cow = 0
foreach t_group in `cohort_groups' {
local cow = `cow'  + 1
qui: sum `beta_hat' if  (`arg_cohort_treatment_date'  == (`t_group')) & (`arg_time' == (`t_group' + `t_ell')) & `_sample'  
if r(N)>0{
qui: replace influence_final = influence_final  + ((r(N)/total_final)*(phi_`cow'_`t_ell') + r(mean)*psi_`cow'_`t_ell')  if  `_sample'
		 }	
	}
}
qui: gen  influence_sqr_final = (influence_final)^2 if `_sample'
qui: sum influence_sqr_final if `_sample'
qui: gen se_theta_final = sqrt(r(mean)/r(N))
/*   REPORT AVERAGE TREATMENT EFFECTS      */
** report to eclass 
ereturn clear 
qui: sum `_sample' if `_sample'
scalar obs_actual = r(N)
qui: sum `copy' if `copy'
scalar N = r(N)
cap : drop sample 
cap : gen sample = `_sample' // actual sample used, made available directly
ereturn post, esample(`copy') obs(`=N') depname(`y')
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
qui: sum `beta_hat' if `_sample' & (`arg_cohort_treatment_date' != 0) & (`arg_time' == (`arg_cohort_treatment_date' + `t_ell' )) 
scalar beta_error = beta_error + r(N)==0
matrix b[`t_ell'+1, 1] = r(mean)
qui: sum se_theta_`t_ell'
matrix se[`t_ell'+1, 1] = r(mean)
}
* Fill with avg. effect (theta)
qui: sum `beta_hat' if `_sample' & (`arg_cohort_treatment_date' != 0) & (`arg_time' >= (`arg_cohort_treatment_date')) 
matrix b[Ncoef, 1] = r(mean)
qui: sum se_theta_final
matrix se[Ncoef, 1] = r(mean)
* Prepare to export 
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
ereturn scalar N_obs = obs_actual
ereturn scalar keep = `keep'
ereturn scalar first = `first'
ereturn scalar periods = `periods'
di "Number of Observations: " obs_actual
end
