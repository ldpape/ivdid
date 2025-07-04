# Implementation in Stata of Instrumented Differences-in-Differences for Staggered Rollout with Heterogeneous Treatment Effects
Before using this package, please consult Miyaji (2025, working paper) at: https://arxiv.org/abs/2405.12083

### Installation of beta version (subject to changes, use at your own risk)

        cap ado uninstall sivdid
        net install sivdid, from("https://raw.githubusercontent.com/ldpape/ivdid/refs/heads/main/")

## Example (see help file)

Load Oreopolous (UK) data :

        use "https://raw.githubusercontent.com/ldpape/ivdid/master/oreopoulos_sample.dta", replace

Generate the cohort treatment date variable by identifying first time period of treatment. Here, there is only one cohort:

        gen tempvar = yearat14*drop15
        replace tempvar = . if tempvar == 0
        egen cohort_treatment_date = min(tempvar)

Identify never treated group:

        replace cohort_treatment_date = 0 if nireland

Example of package usage :

        sivdid, y(learn) d(agelfted) z(drop15) cohort_treatment_date(cohort_treatment_date) time(yearat14) periods(9) permanent keep(2)

Generate the event study plot and use with "if" option (blue line is the estimated weighted average) :

        sivdid if learn>0 , y(learn) d(agelfted) z(drop15) cohort_treatment_date(cohort_treatment_date) time(yearat14) periods(9) permanent keep(2) event_study
        
   ![oreo_graph](https://github.com/user-attachments/assets/3069fdd6-dd49-44d9-ae01-5024f62894f3)

        graph export "ivdid_eventstudy.pdf",replace 

### Export estimates using estout 

           eststo: sivdid, y(learn) d(agelfted) z(drop15) cohort_treatment_date(cohort_treatment_date) time(yearat14)  periods(9)
           esttab *
### Bootstrapping for endogenous count or control variables
Use with endogenous count model of Mullahy (1996) : 

        gen earn = exp(learn) 
        sivdid, y(earn) d(agelfted) z(drop15) cohort_treatment_date(cohort_treatment_date) time(yearat14) periods(9) permanent keep(2)  exponential event_study
            graph_options(xtitle("New x-axis title")
        eststo: bootstrap _b ,  level(95)   seed(1234) reps(10)  : sivdid, y(earn) d(agelfted) z(drop15) cohort_treatment_date(cohort_treatment_date) time(yearat14)  periods(9) exponential

### Citations
Please cite the original contribution by Sho Miyaji (2025) and applied paper for which this package was originally developped:

@misc{miyaji2025instrumenteddifferenceindifferencesheterogeneoustreatment,
      title={Instrumented Difference-in-Differences with Heterogeneous Treatment Effects}, 
      author={Sho Miyaji},
      year={2025},
      eprint={2405.12083},
      archivePrefix={arXiv},
      primaryClass={econ.EM},
      url={https://arxiv.org/abs/2405.12083}, 
}

@techreport{helmers2025judge,
  author       = {Helmers, Christian and Love, Brian J. and Pape, Louis-Daniel},
  title        = {Judge (Ideology) Shopping},
  institution  = {Santa Clara Univ. Legal Studies Research Paper No. 5143777},
  year         = {2025},
  month        = {February},
  note         = {Available at SSRN: \url{https://ssrn.com/abstract=5143777}}
}
