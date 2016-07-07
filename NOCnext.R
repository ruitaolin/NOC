##########################################################################
#          Implementation of  Nonparametric Overdose Control             #
#                    Lin Ruitao and Yin Guosheng                         #
#                            gyin@hku.hk                                 #
##########################################################################
# -----------------------------------------------------------------------#
#     Main calculation function                                          #
#     get.next.noc: generate next dose level for conventional cases      #
#                                                                        #
#                             ---Input---                                #
#             target --> target toxicity probability                     #
#             dlt -->  toxicity outcome for each patient                 #
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
#            * A summary table                                           #
# -----------------------------------------------------------------------#
##########################################################################

 
get.next.noc <- function(target, dlt, dose.level,
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
        y = rep(0, ndose)
        n = rep(0, ndose)
		    for (j in 1:ndose)
		    {
		     y[j] = sum(dlt[which(dose.level==j)])
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
		
		
		pos.model<-round(pos.model,2)
		cat("The feasibility bound for overdose control rule is", a, "\n")
		cat("The dose-switching cutoff is", eta, "\n")
		cat("The dose-elimination cutoff is", lambda, "\n")
		cat("The posterior model probabilities are", pos.model, "\n")
		cat("The posterior probability that the current dose level is overly toxic is", pos.tox, "\n")
		cat("The next dose level is", d, "\n")
		summary = cbind(seq(1,length(dose.level)), dose.level, dlt)
		colnames(summary) = c("Patient No.", "Dose level", "DLT")
		return(summary)
}


# generate the next dose level after the 4th in the sonidegib trial
# by treating the missing data as no toxicity.
target = 0.33
ndose = 5
dlt= c(0,0,0,0,0,0,1,0,0,1,0,0)
dose.level=c(1,1,1,2,2,2,3,3,3,3,3,3)

get.next.noc(target, dlt, dose.level, ndose)


