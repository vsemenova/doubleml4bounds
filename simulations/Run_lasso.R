library(hdm)
library(xtable)

directoryname<-"/n/tata/doubleml4bounds/"
setwd(directoryname)

source("simulations/Functions.R")
source("simulations/FirstStage.R")
source("simulations/RemoveControls.R")


## true parameter
beta = c(1,1)
## true Sigma
Sigma = matrix (c(1,0.5,0.5,1),nrow=2,ncol=2)
## true regression parameters
decay=2
cons.x=c(2,1)
cons.y=1
cons.eta=1
rho=0.5
p=50


## number of points on the boundary
num_circle=50
N_reps=500
num_splits<-2
B=500
control.name = paste0("Z.",as.character(1:p))
### confidence level
ci_alpha=0.05


## save results here
sample_sizes=c(250,500,700,1000)
Deltas = c(1,2,3)

total_risk<-array(0,c(2,length(sample_sizes),length(Deltas)))
outer_risk<-array(0,c(2,length(sample_sizes),length(Deltas)))
inner_risk<-array(0,c(2,length(sample_sizes),length(Deltas)))
rej.freq<-array(0,c(2,length(sample_sizes),length(Deltas)))

### save results here
finaltable<-matrix(0,nrow=length(sample_sizes),ncol=length(Deltas)*4)

### first stage estimation
method.treat=rlasso
method.outcome=rlasso

### critical value based on bootstrap
for (j in 1:length(Deltas)) {
  ### true value of support function is approximated by simulation
  mydata<-simulate_data_2d(seed=1, decay=decay,
                           sample_size=100000,p=p,beta=beta,cons.x=cons.x,
                           cons.eta=cons.eta,
                           Sigma=Sigma,Delta=Deltas[j],rho=rho,control.name=control.name)
  
  
  res<-estimate_sigma(num_circle=num_circle,Delta=Deltas[j],treat.name=c("V1","V2"),
                      outcome.name=c("yL"),is_list=FALSE,mydata, sample_size=100000,partial_out=FALSE)
  for (i in 1:length(sample_sizes)) {
    print(paste0("Sample size is ",sample_sizes[i]))
    print(paste0("Bracket width is ", Deltas[j]))
    
    
    sample_size=sample_sizes[i]
    #set.seed(888)
  #  indices<-sample(1:sample_size,size=sample_size,replace=FALSE)
    #INDS<-lapply(1:num_splits, function (i) indices[(ceiling(sample_size/num_splits)*(i-1) + 1): min((ceiling(sample_size/num_splits)*(i)),sample_size)])
    
    
    simulated_data<-lapply(1:N_reps,simulate_data_2d, decay=decay,
                           sample_size=sample_sizes[i],p=p,beta=beta,cons.x=cons.x,
                           cons.eta=cons.eta,
                           Sigma=Sigma,Delta=Deltas[j],rho=rho,control.name=control.name)
    
    fitted_residuals<-lapply(simulated_data,first_stage,method.treat =method.treat,
                             method.outcome = method.outcome,
                             treat.name = c("X1","X2"),
                             outcome.name = c("yL"),
                             control.name = control.name,num_splits=1)
    
    
    
    
    estimated_support_function<-lapply(fitted_residuals,estimate_sigma,sample_size=sample_sizes[i],num_circle=num_circle,
                                       Delta=Deltas[j],treat.name=c("treat.X1","treat.X2"),outcome.name=c("true_outcome"),fitted_outcome_name="fitted_outcome",
                                       partial_out=TRUE,weighted=FALSE)
    
    bootstrap_support_function<-lapply(1:B,estimate_sigma,estimate_sigma, sample_size=sample_sizes[i],num_circle=num_circle,
                                       Delta=Deltas[j],treat.name=c("treat.X1","treat.X2"),outcome.name=c("true_outcome"),fitted_outcome_name="fitted_outcome",
                                       partial_out=TRUE,data=fitted_residuals[[1]],weighted=FALSE)
    
    diff<-array(0,c(num_circle,2,N_reps))
    for (k in 1:N_reps) {
      diff [,,k] = estimated_support_function[[k]] - res
      
    }
    ## maximum estimated error over the boundary
    sup_stat=apply(diff,c(2,3),function(x) max(abs(x)))
    
    
    
    ## maximum mistake over the boundary
    total_risk[,i,j]<-apply(sup_stat,1,mean )
    outer_risk[,i,j]<-apply(apply(diff,c(2,3),function(x) max(vabs(x))),1,mean )
    inner_risk[,i,j]<-apply(apply(diff,c(2,3),function(x) max(vabs(-x))),1,mean )
    
    ## bootstrap critical value
    diff_boot<-array(0,c(num_circle,2,B))
    for (b in 1:B) {
      diff_boot [,,b] = bootstrap_support_function[[b]] - estimated_support_function[[1]]
      
    }
    sup_stat_boot=apply(diff_boot,c(2,3),function(x) max(abs(x)))
    crit.value<-apply(sup_stat_boot,1,quantile,1-ci_alpha)
    
    ## rejection frequency
    rej.freq[,i,j]<-apply(sup_stat>crit.value,1,mean )
    
    finaltable[i,1:4+4*(j-1)] <-c(total_risk[1,i,j],outer_risk[1,i,j],inner_risk[1,i,j],rej.freq[1,i,j])
    
  }
  
  
}

finaltable<-apply(finaltable,2,round,2)
write.csv(finaltable,paste0(directoryname,"/simulations/Tables/Table_lasso1.csv"))
