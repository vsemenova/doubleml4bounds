library(hdm)
library(xtable)

directoryname<-"/Users/virasemenova/Dropbox (MIT)/JMP_Submission/JoE/code/"
setwd(directoryname)

source("empirical_application/Functions.R")
mydata<-read.csv("data/cps2015.csv")
mydata<-mydata[1:10000,]
### save results here

res_table<- matrix(0, 6, 2)	

#### POINT-IDENTIFIED SPECIICATION ####


Z = model.matrix(~ -1 + (exp1+exp2+exp3+exp4)*(educ+occ2+ind2+ms+reg)
                 +(educ+occ2+ind2+ms+reg),data=mydata )
dim(Z)

center_apply <- function(x) {
  apply(x, 2, function(y) y - mean(y))
}


Z = center_apply(Z) #'recenter contols; this is mostly not needed


library(hdm)

fmla.y<- mydata$lnw~ Z
fmla.d<- mydata$female~ Z

rY<- rlasso(fmla.y)$res	
rD<- rlasso(fmla.d)$res	
partial.fit.lasso<- lm(rY~rD-1)	
HCV.coefs <- vcovHC(partial.fit.lasso, type = 'HC');
partial.lasso.se <- sqrt(diag(HCV.coefs))
res_table[1,1]<-coef(partial.fit.lasso)
res_table[1,2]<-partial.lasso.se
#### MID-POINT SPECIICATION ####



####

Deltas<-c(1,1)

for (Delta in Deltas) {
      
   ### create bracketed wage 
  bracketed_wage<-bracket(mydata$lnw,Delta=Delta,1,6)
  
  bracketed_wage$yM<-(bracketed_wage$yL+bracketed_wage$yU)/2
  rY.M<- rlasso(x=Z,y=bracketed_wage$yM)$res
  partial.fit.lasso.MID<- lm(rY.M~rD-1)	
  HCV.coefs <- vcovHC(partial.fit.lasso, type = 'HC');
  partial.lasso.se <- sqrt(diag(HCV.coefs))
  
  res_table[2,1]<-coef(partial.fit.lasso.MID)
  res_table[2,2]<-partial.lasso.se
  
  
  
  
  ### bounds specification (the simplest possible)
  
  
  
  covariate = rD
  YL.fit<- rlasso(x=Z,y=bracketed_wage$yL)
  rYL<-YL.fit$res
  
  outcomeU = rYL + Delta*(rD>0) - Delta/2
  partial.fit.lasso.U<- lm(outcomeU~rD-1)	
  HCV.coefs.U <- vcovHC(partial.fit.lasso.U, type = 'HC');
  partial.lasso.se.U <- sqrt(diag(HCV.coefs.U)); #'White std errors
  res_table[3,1]<-coef(partial.fit.lasso.U)
  res_table[3,2]<-partial.lasso.se
    
  outcomeL = rYL + Delta*(rD<0) - Delta/2
  partial.fit.lasso.L<- lm(outcomeL~rD-1)	
  HCV.coefs.L <- vcovHC(partial.fit.lasso.L, type = 'HC');
  partial.lasso.se.L <- sqrt(diag(HCV.coefs.L)); #'White std errors
  res_table[4,1]<-coef(partial.fit.lasso.L)
  res_table[4,2]<-partial.lasso.se.L 
  
  ### bounds specification (the hardest possible)
  
  covariate = rD
  
  Sign<-rlassologit(x=Z,y=(covariate>0))
  rSign<-predict(Sign, data=Z,type="response")
  outcomeU.2 = rYL + Delta*(rD>0) - Delta*rSign
  partial.fit.lasso.U.2<- lm(outcomeU.2~rD-1)	
  HCV.coefs.U.2 <- vcovHC(partial.fit.lasso.U.2, type = 'HC');
  partial.lasso.se.U.2 <- sqrt(diag(HCV.coefs.U.2)); #'White std errors
  res_table[5,1]<-coef(partial.fit.lasso.U.2)
  res_table[5,2]<-partial.lasso.se.U.2
  
  
  outcomeL.2 = rYL + Delta*(rD<0) - Delta*(1-rSign)
  partial.fit.lasso.L.2<- lm(outcomeL.2~rD-1)	
  HCV.coefs.L.2 <- vcovHC(partial.fit.lasso.L.2, type = 'HC');
  partial.lasso.se.L.2 <- sqrt(diag(HCV.coefs.L.2))
  res_table[6,1]<-  coef(partial.fit.lasso.L.2)
  res_table[6,2]<- partial.lasso.se.L.2
}

colnames(table)<- names(summary(full.fit)$coef["female",])[1:2]	
rownames(table)<- c("without controls","full reg", "partial reg")	
tab<- xtable(table, digits=c(2, 5,10))	



