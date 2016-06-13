##########################################################################
#          Implementation of  Nonparametric Overdose Control             #
#                    Lin Ruitao and Yin Guosheng                         #
#                            gyin@hku.hk                                 #
##########################################################################
# -----------------------------------------------------------------------#
#     Main calculation function                                          #
#     select.mtd.noc: select the MTD based on the full data              #
#                                                                        #
#                             ---Input---                                #
#             target --> target toxicity probability                     #
#             dlt -->  toxicity outcome for each patient                 #
#             dose.level --> dose level for each patient                 #
#             ndose --> number of dose levels                            #
#             epi --> a small number that defines the neighbour of target#
#             lambda --> dose-elimination cutoff                         #
#                                                                        #
#                             ---Output---                               #
#            * posterior model probabilities                             #
#            * MTD recommendation                                        #
# -----------------------------------------------------------------------#
##########################################################################

 
select.mtd.noc <- function(target, dlt, dose.level,
 ndose, epi=0.05, lambda=0.85){ 
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
			pos.model[elimi==1]=0;
			d_opt<-which.max(pos.model)

		
		pos.model<-round(pos.model,2)
		cat("The posterior model probabilities are", pos.model, "\n")
		cat("The MTD is the dose level", d_opt, "\n")

}


# generate the next dose level after the 4th in the sonidegib trial
target = 0.33
ndose = 5
dlt= c(0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,1,1,1,0,1,0,0,0,1,0)
dose.level=c(1,1,1,2,2,2,3,3,3,3,3,3,2,2,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2)

select.mtd.noc(target, dlt, dose.level, ndose)


