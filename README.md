# NOC (Nonparametric Overdose Control)
R codes to implement nonparametric overdose control and fractional nonparametric overdose control designs in phase I dose-finding trials.
# Description
The primary objective of phase I oncology trials is to identify the maximum tolerated dose (MTD) whose induced dose-limiting toxicity (DLT) probability is closest to the target toxicity rate. Under the framework of Bayesian model selection, we propose a nonparametric overdose control (NOC) design for dose finding in phase I clinical trials. Each dose assignment is guided via a feasibility bound, which thereby can control the number of patients allocated to excessively toxic dose levels. Operatively, the NOC design requires that the DLT outcome is ascertainable quickly after treatment. However, many types of DLTs do not exhibit till a long period of time after treatment in real applications. As a result, the dose-escalation procedure of the trial would be affected due to such late-onset toxicities. We further propose a fractional NOC (fNOC) design in conjunction with a so-called fractional imputation approach, to account for late-onset toxicity outcomes.
#Functions
The repository includes three functions:
* NOCnext.R: The R code that includes the function ```get.next.noc``` to select the next dose level for the new patients by the NOC design.
```rscript
get.next.noc(target, dlt, dose.level, ndose, epi, a, eta, lambda)
```
* fNOCnext.R: The R code that include the function ```get.next.fnoc``` to select the next dose level for the new patents by the fNOC design when the toxicity outcome is late-onset.
```rscipt
get.next.fnoc(target, enter.time, dlt.time, current.time, tau, dose.level, ndose, epi, a, eta, lambda)
```
* NOCmtd.R: The R code to estimate the MTD level at the end of the trial.
* 

#Inputs
* ```target```: The target toxicity probability, e.g., ```target=0.33```.
* 

#Authors and Reference
* Ruitao Lin and Guosheng Yin (gyin@hku.hk)
* Lin, R. and Yin, G. (2016) “Nonparametric overdose control with late-onset toxicity in phase I clinical trials”

