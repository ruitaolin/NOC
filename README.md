# NOC (Nonparametric Overdose Control)
R codes to implement nonparametric overdose control and fractional nonparametric overdose control designs in phase I dose-finding trials.
# Description
The primary objective of phase I oncology trials is to identify the maximum tolerated dose (MTD), whose induced dose-limiting toxicity (DLT) probability is the closest to the target toxicity rate. Under the framework of Bayesian model selection, we propose a nonparametric overdose control (NOC) design for dose finding in phase I clinical trials. Each dose assignment is guided via a feasibility bound, which thereby can control the number of patients allocated to excessively toxic dose levels. Operatively, the NOC design requires that the DLT outcome is ascertainable quickly after treatment. However, many types of DLTs do not exhibit till a long period of time after treatment in real applications. As a result, the dose-escalation procedure of the trial would be affected due to such late-onset toxicities. We further propose a fractional NOC (fNOC) design in conjunction with a so-called fractional imputation approach, to account for late-onset toxicity outcomes.
#Functions
The repository includes three functions:
* NOCnext.R: The R code that includes the function ```get.next.noc``` to select the next dose level for the new patients by the NOC design.
```rscript
get.next.noc(target, dlt, dose.level, ndose, epi, a, eta, lambda)
```
* fNOCnext.R: The R code that includes the function ```get.next.fnoc``` to select the next dose level for the new patents by the fNOC design when the toxicity outcome is late-onset.
```rscipt
get.next.fnoc(target, enter.time, dlt.time, current.time, tau, dose.level, ndose, epi, a, eta, lambda)
```
* NOCmtd.R: The R code that includes the function ```select.mtd.noc``` to estimate the MTD level at the end of the trial.
```rscript
select.mtd.noc(target, dlt, dose.level, ndose, epi, lambda)
```


#Inputs
* ```target```: The target toxicity probability, e.g., ```target<-0.33```.
* ```dlt```: A vector of length *n* that stores the toxicity outcome for each patient, where *n* is the total number of patients so far.
* ```dose.level```: A vector of length *n* that stores the dose level assigned to each patient.
* ```ndose```: Number of prespecified dose levels of the new drug.
* ```epi```: A small positive value that defines the neighbourhood of the target toxicity probability.
* ```a```: The feasibility bound for overdose control, as default, ```a<-0.35```. 
* ```eta```: The dose-switching cutoff, as default, ```eta<-0.60```.
* ```lambda```: The dose-elimination cutoff, as default, ```lambda<-0.85```.
* ```enter.time```: A vector of length *n* that stores the day of arrival of each patient under late-onset cases.
* ```dlt.time```: A vector of length *n* that stores the time-to-toxicity outcome of each patient; If the subject has not experienced the DLT by the decision-making time, his time-to-toxcity outcome is *0*.
* ```current.time```: The arrival time of the new patient, or the decision-making time to decide the next dose level.
* ```tau```: The length of follow-up period.


#Example
We apply the NOC and fNOC designs to the sonidegib trial.
* Based on the accumulated data, two DLTs were observed at dose level 3 when patient 13 arrived on day 130. At this moment, patients 6, 8, 9, 11, and 12 were still under the follow-up of evaluation without experiencing any DLT, which led to a total of five missing toxicity outcomes. We utilize the following code to decide the dose level for patient 13.
```rscript
target <- 0.33
ndose <- 5
enter.time <- c(4,7,19,29,31,50,58,67,78,91,100,118)
dlt.time <- c(0,0,0,0,0,0,65,0,0,29,0,0)
dose.level <- c(1,1,1,2,2,2,3,3,3,3,3,3)
tau <- 90
current.time <- 130
get.next.fnoc(target, enter.time, dlt.time, current.time,tau, dose.level, ndose)
```
The output is given by 
```rscript
The feasibility bound for overdose control rule is 0.35 
The dose-switching cutoff is 0.6 
The dose-elimination cutoff is 0.85 
The posterior model probabilities are 0.02 0.16 0.55 0.2 0.07 
The posterior probability that the current dose level is overly toxic is 0.4749329 
The next dose level is 2 
      Patient No. Dose level Day of arrival Observed DLT Time to DLT Fractional DLT
 [1,]           1          1              4            0          91      0.0000000
 [2,]           2          1              7            0          91      0.0000000
 [3,]           3          1             19            0          91      0.0000000
 [4,]           4          2             29            0          91      0.0000000
 [5,]           5          2             31            0          91      0.0000000
 [6,]           6          2             50            0          91      0.0000000
 [7,]           7          3             58            1          65      1.0000000
 [8,]           8          3             67            0          91      0.1428571
 [9,]           9          3             78            0          91      0.1428571
[10,]          10          3             91            1          29      1.0000000
[11,]          11          3            100            0          91      0.1428571
[12,]          12          3            118            0          91      0.2207792
```
* On the other hand, if we had treated these missing data as no DLTs, we can utilize the ```get.next.noc``` function to generate the next dose level, as given by
```rscript 
target <- 0.33
ndose <- 5
dlt <- c(0,0,0,0,0,0,1,0,0,1,0,0)
dose.level <- c(1,1,1,2,2,2,3,3,3,3,3,3)
get.next.noc(target, dlt, dose.level, ndose)

-----------------------output------------------------
The feasibility bound for overdose control rule is 0.35 
The dose-switching cutoff is 0.6 
The dose-elimination cutoff is 0.85 
The posterior model probabilities are 0.01 0.08 0.49 0.28 0.14 
The posterior probability that the current dose level is overly toxic is 0.3341323 
The next dose level is 3 
      Patient No. Dose level DLT
 [1,]           1          1   0
 [2,]           2          1   0
 [3,]           3          1   0
 [4,]           4          2   0
 [5,]           5          2   0
 [6,]           6          2   0
 [7,]           7          3   1
 [8,]           8          3   0
 [9,]           9          3   0
[10,]          10          3   1
[11,]          11          3   0
[12,]          12          3   0
```
* In the end, three dose levels had been explored with four DLTs out of 18 patients at dose level 2 and four DLTs out of 9 patients at dose level 3 being observed. The MTD can be estimated by the ```select.mtd.noc``` function. 
```rscript 
target <- 0.33
ndose <- 5
dlt <- c(0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,1,1,1,0,1,0,0,0,1,0)
dose.level <- c(1,1,1,2,2,2,3,3,3,3,3,3,2,2,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2)
select.mtd.noc(target, dlt, dose.level, ndose)

-----------------------output------------------------
The posterior model probabilities are 0.03 0.55 0.36 0.05 0.01 
The MTD is the dose level 2 
```
# Authors and Reference
* Ruitao Lin and Guosheng Yin
* Lin, R. and Yin, G. (2017) “Nonparametric overdose control with late-onset toxicity in phase I clinical trials”, Biostatistics, 18: 180-194.

