{smcl}
{* *! version 1.0 27june2025}{...}
{viewerjumpto "Description" "sivdid##description"}{...}
{viewerjumpto "Syntax" "sivdid##syntax"}{...}
{viewerjumpto "Citation" "sivdid##citation"}{...}
{viewerjumpto "Authors" "sivdid##contact"}{...}
{viewerjumpto "Examples" "sivdid##examples"}{...}
{viewerjumpto "Stored results" "sivdid##results"}{...}

{title:Title}

{p2colset 5 18 20 2}
{p2col :{cmd:sivdid} {hline 2}} Instrumented Differences-in-differences with Heterogeneous Treatment Effects {p_end}
{p2colreset}

{marker description}{...}
{title:Description}

{pstd} {cmd:sivdid} implements instrumented differences-in-differences with heterogeneous treatment effects and staggered roll-out, as described by {browse "https://arxiv.org/abs/2405.12083": Sho Miyaji (2025, working paper)}. The package first estimates the local average treatment effets (LATE) for all two time periods and two groups (2x2) cases. It then reports a weighted average of these treatment effects where the weights depend on the share of treated in each cohort.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:sivdid}  {ifin}, {cmd:y(outcome variable)}   {cmd:d(endogenous variable)}  {cmd:z(instrumental variable)}  {cmd:cohort_treatment_date(cohort identifier)}  {cmd:time(time identifier)}  {cmd:periods(number of treatment periods)}   {cmd:keep(integer)}  {cmd:permanent}   {cmd:controls(control  variables)}  {cmd:event_study}  {cmd:graph_options(tw options)}  {cmd:exponential}  {cmd:first(it:real number)}   [{help sivdid##options:options}] {p_end}



{synoptset 22}{...}
{synopthdr:Variables}
{synoptline}
{synopt:{it:y}} Outcome variable. {p_end}
{synopt:{it:d}} Endogenous variable. {p_end}
{synopt:{it:z}} Instrumental variable. {p_end}
{synopt:{it:cohort_treatment_date}} Indicates the time period in which individual gets treated. Your chosen control group should have a zero value. {p_end}
{synopt:{it:time}} Time variable. {p_end}
{synopt:{it:periods}} Number of post-treatment periods. {p_end}
{synoptline}

{marker opt_summary}{...}
{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt permanent}{cmd:(}{help sivdid##permanent:permanent}{cmd:)}} Drop cohort which are not present the full time horizon (defined by periods). {p_end}

{synopt:{opt keep}{cmd:(}{help sivdid##keep:keep}{cmd:)}} Drop cohort who don't have at least this number of observations at each moment in time. Default set at 0. Requires the permanent option. {p_end}

{synopt:{opt controls}{cmd:(}{help sivdid##controls:controls}{cmd:)}} You can add a varlist of control variables. Note that the standard errors do not account for them. It may be a good idea to bootstrap your standard errors in this case. {p_end}

{synopt:{opt event_study}{cmd:(}{help sivdid##event_study:event_study}{cmd:)}} Create the event-study graph for the per period average LATE. The blue horizontal line is the estimated total average. The graph includes 95% confidence intervals. {p_end}

{synopt:{opt graph_options}{cmd:(}{help sivdid##graph_options:graph_options}{cmd:)}} You can write here options that are then provided to two way graphs. {p_end}

{synopt:{opt exponential}{cmd:(}{help sivdid##exponential:exponential}{cmd:)}} Point estimates calculated based on using an endogenous count model. Standard errors are not provided and need to be bootstrapped. {p_end}

{synopt:{opt first}{cmd:(}{help sivdid##exponential:exponential}{cmd:)}} Not commonly used but allows you to remove 2x2 DiD cases where the F-test is below the value you set here. Note that this creates selection bias which invalidates inference based on the provided standard errors. {p_end}

{marker postestimation}{...}
{title:Post-Estimation}

{pstd} The program generates two variables:{p_end}
{phang2} {cmd:_late_hat} is the estimated LATEs at the 2x2 DiD level with associated standard errors in  {cmd:_late_hat_se} 

{marker citation}{...}
{title:Citation}

{pstd}

This package was developed for the following research: 
C. Helmers, B. Love and L-D Pape, Judge (Ideology) Shopping, {browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5143777": Santa Clara Univ. Legal Studies Research Paper No. 5143777}

This package is based on the following article: 
Sho Miyaji, Instrumented Difference-in-Differences with Heterogeneous Treatment Effects (2025). 
{browse "https://arxiv.org/abs/2405.12083":arXiv, 2405.12083}. 

{pstd}BibTeX:{p_end}
@misc{miyaji2025,
      title={Instrumented Difference-in-Differences with Heterogeneous Treatment Effects}, 
      author={Sho Miyaji},
      year={2025},
      eprint={2405.12083},
      archivePrefix={arXiv},
      primaryClass={econ.EM},
      url={https://arxiv.org/abs/2405.12083}, 
}

{marker examples}{...}
{title:Examples (Oreopolous Data)}

Load Oreopolous (UK) data :
{phang2}{cmd: use "https://raw.githubusercontent.com/ldpape/ivdid/master/oreopoulos_sample.dta", replace}{p_end}

Generate the cohort treatment date variable by identifying first time period of treatment.
Here, there is only one cohort:

{phang2}{cmd: gen tempvar = yearat14*drop15 }{p_end}
{phang2}{cmd: replace tempvar = . if tempvar == 0 }{p_end}
{phang2}{cmd: egen cohort_treatment_date = min(tempvar) }{p_end}

Identify never treated group:
{phang2}{cmd:  replace cohort_treatment_date = 0 if nireland }{p_end}

Example of package usage :
{phang2}{cmd:   sivdid, y(learn) d(agelfted) z(drop15) cohort_treatment_date(cohort_treatment_date) time(yearat14)  periods(9)   }{p_end}
{phang2}{cmd:   sivdid, y(learn) d(agelfted) z(drop15) cohort_treatment_date(cohort_treatment_date) time(yearat14)  periods(9) permament keep(2)   }{p_end}
{phang2}{cmd:   sivdid, y(learn) d(agelfted) z(drop15) cohort_treatment_date(cohort_treatment_date) time(yearat14)  periods(9) permament keep(2)  event_study  }{p_end}
{phang2}{cmd:   sivdid if learn>0 , y(learn) d(agelfted) z(drop15) cohort_treatment_date(cohort_treatment_date) time(yearat14)  periods(9) permament keep(2)  event_study graph_options(xtitle("New x-axis title")  }{p_end}

{marker results}{...}
{title:Stored Results}

{pstd}{cmd:sivdid} stores the following in {cmd:e()}:{p_end}

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}} Number of observations {p_end}
{synopt:{cmd:e(sample)}} Sample used for estimation {p_end}
{synopt:{cmd:e(late)}} Matrix of per period average LATEs {p_end}
{synopt:{cmd:e(std)}} Matrix of standard errors of per period average LATEs {p_end}
{synopt:{cmd:e(error)}} Equal to one if an average per period LATE could not be identified {p_end}
{synopt:{cmd:e(periods)}}  Number of selected periods {p_end}
{synopt:{cmd:e(first)}} F-test restriction (seldom used) {p_end}
{synopt:{cmd:e(keep)}} Indicated permanent cohort restriction {p_end}

{marker author}{...}
{title:Author}

{pstd} Louis Pape {break}
Télécom Paris (CREST - IP Paris) {break}
Contact: {browse "mailto:louis.pape@telecom-paris.fr":louis.pape@telecom-paris.fr}{p_end}
