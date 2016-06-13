##########################################################################
#    Implementation of Fractional Nonparametric Overdose Control         #
#                    Lin Ruitao and Yin Guosheng                         #
#                            gyin@hku.hk                                 #
##########################################################################
# -----------------------------------------------------------------------#
#     Main calculation function                                          #
#     get.next.fnoc: generate next dose level for late-onset cases       #
#                                                                        #
#                             ---Input---                                #
#             target --> target toxicity probability	                 #
#             enter.time --> days of arrival for each patient            #
#             dlt.time  -->  time to dlt for each patient                #
#                            0 indicates no dlt                          #
#             current.time --> day of arrival for the the patient        #
#             tau --> length of follow-up period                         #
#             dose.level --> dose level for each patient                 #
#             ndose --> number of dose levels                            #
#             epi --> a small number that defines the neighbour of target#
#             a --> feasibility bound for overdose control               #
#             eta --> dose-switching cutoff                              #
#             lambda --> dose-elimination cutoff                         #
#                                                                        #
#                             ---Output---                               #
#            * posterior model probabilities                             #
#            * posterior probability that the current dose is too toxic  #
#            * next dose recommendation                                  #
#            * Plot of Kaplan-Meier estimate                             #
#            * A summary table                                           #
# -----------------------------------------------------------------------#
##########################################################################

library(survival)

get.next.fnoc <- function(target, enter.time, dlt.time, current.time, tau, dose.level,
 ndose, epi=0.05, a=0.35, eta=0.6, lambda=0.85){ 
posteriorH<-function(y,n,target,p.sample,d){
	lik<-rep(1,length(y))
	pos.tox<-rep(0,length(y))
	for (i in 1:length(y)){
		likeli<-1;
			for (j in 1:ndose){
				likeli<-likeli*p.sample[,j,i]^y[j]*(1-p.sample[,j,i])^(n[j]-y[j])
			}

		lik[i]<-mean(likeli)
		pos.tox[i]<-sum(likeli[p.sample[,d,i]>target])/sum(likeli)
	}
lik<-lik/sum(lik)
pos.tox<-sum(pos.tox*lik)	
return(c(lik,pos.tox))
}

set.seed(6)	
NN<-50000
p.sample<-array(matrix(nrow=10,ncol=6),dim=c(NN,ndose,ndose))
for (i in 1:ndose){
	p.sample[,i,i]<-runif(NN,target-epi,target+epi)
	if (i>1){
		for (j in (i-1):1){
			 p.sample[,j,i]<-runif(NN,0,sapply(p.sample[,j+1,i],function(x) min(x,target-epi	)))
		}
	}
	if (i<ndose){
		for (j in (i+1):ndose){
			p.sample[,j,i]<-runif(NN,sapply(p.sample[,j-1,i],function(x) max(x,target+epi)),0.8)
		}
	}
	
}
		d=dose.level[length(dose.level)] 
		elimi = rep(0, ndose); 		
        assesstime = enter.time+tau;	
        dlt.time[dlt.time==0]= tau+1;
        yo = (dlt.time<=tau)*(assesstime<=current.time)+(dlt.time<=(current.time+tau-assesstime))*(current.time<assesstime);		
		if (sum(yo)==0)	{stop("Must observe at least one DLT")}
	    if (sum(yo)!=0){			
			
			otime = yo*dlt.time+(1-yo)*((current.time+tau-assesstime)*(current.time<assesstime)+tau*(assesstime<=current.time))			
			kmfit = survfit(Surv(otime,yo)~1)	
			plot(kmfit,xlab="Time (days)",ylab="Survival probability",
			main="Kaplan-Meier estimate",cex.lab=1.3,cex.main=1.3,cex.sub=1.3)	
			ym = yo
			
			for (i in 1:length(yo)){
			if (current.time<assesstime[i] & yo[i]==0){
			ym[i]=(kmfit$surv[tail(which(kmfit$time<=(current.time-assesstime[i]+tau+0.001)),n=1)]-
			kmfit$surv[tail(which(kmfit$time<=tau),n=1)])/
			kmfit$surv[tail(which(kmfit$time<=(current.time-assesstime[i]+tau+0.001)),n=1)]
			}
			}
            
            y=rep(0,ndose)
			n=y;
		    for (j in 1:ndose)
		    {
		     y[j] = sum(ym[which(dose.level==j)])
		     n[j] = sum(dose.level==j)
		    }

			pos.model=posteriorH(y,n,target,p.sample,d)
			pos.tox<-pos.model[ndose+1];
			if(pos.tox>lambda)
		    ## determine if the current dose should be eliminated
			{
					elimi[d:ndose]=1;
					if (elimi[1]==1) {stop("The trial should be terminated early"); }
			}			
			pos.model<-pos.model[1:ndose];

			pos.model=pos.model/sum(pos.model)
			d_opt<-which.max(pos.model)
			if(pos.model[d_opt]<eta){
			pos.cum.model<-cumsum(pos.model)
            d_opt<-which.min(abs(pos.cum.model-a))
            }
 
		 if(d < d_opt && d!=ndose) { if(elimi[d+1]==0) d=d+1; }
			else if(d > d_opt && d!=1) { d=d-1; }
			else { d=d; } 
		
		}
		pos.model<-round(pos.model,2)
		cat("The feasibility bound for overdose control rule is", a, "\n")
		cat("The dose-switching cutoff is", eta, "\n")
		cat("The dose-elimination cutoff is", lambda, "\n")
		cat("The posterior model probabilities are", pos.model, "\n")
		cat("The posterior probability that the current dose level is overly toxic is", pos.tox, "\n")
		cat("The next dose level is", d, "\n")
		summary = cbind(seq(1,length(dose.level)), dose.level, enter.time, dlt.time<=tau, dlt.time, ym)
		colnames(summary) = c("Patient No.", "Dose level", "Day of arrival", "Observed DLT", "Time to DLT", "Fractional DLT")
		return(summary)
}


# generate the next dose level after the 4th in the sonidegib trial
target = 0.33
ndose = 5
enter.time = c(4,7,19,29,31,50,58,67,78,91,100,118)
dlt.time = c(0,0,0,0,0,0,65,0,0,29,0,0)
dose.level=c(1,1,1,2,2,2,3,3,3,3,3,3)
tau = 90
current.time=130
get.next.fnoc(target, enter.time, dlt.time, current.time,tau, dose.level, ndose)


