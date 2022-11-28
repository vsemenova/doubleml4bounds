#install.packages("hdm")
#install.packages("randomForest")
#install.packages("sandwich")
library(hdm)
library(randomForest)
library(sandwich)
library(ranger)
library(xtable)
res_table<- matrix(0, 4, 4)	

directoryname<-"/n/tata/doubleml4bounds/"
setwd(directoryname)

source("empirical_application/Functions.R")
mydata<-read.csv("data/cps2015.csv")
mydata[,c("occ","occ2","ind","ind2","educ","ms","reg","other")]<-sapply(mydata[,c("occ","occ2","ind","ind2","educ","ms","reg","other")],function(x) as.character(x))
mydata[,c("occ","occ2","ind","ind2","educ","ms","reg","other")]<-sapply(mydata[,c("occ","occ2","ind","ind2","educ","ms","reg","other")],function(x) as.factor(x))
mydata$occ<-as.factor(as.character(mydata$occ))
mydata$occ2<-as.factor(as.character(mydata$occ2))
mydata$ind<-as.factor(as.character(mydata$ind))
mydata$ind2<-as.factor(as.character(mydata$ind2))
mydata$educ<-as.factor(as.character(mydata$educ))
mydata$ms<-as.factor(as.character(mydata$ms))
mydata$reg<-as.factor(as.character(mydata$reg))
mydata$other<-as.factor(as.character(mydata$other))

Z = model.matrix(~ -1 + (exp1+exp2+exp3+exp4)*(educ+occ2+ind2+ms+reg)
                 +(educ+occ2+ind2+ms+reg),data=mydata )
dim(Z)
center_apply <- function(x) {
  apply(x, 2, function(y) y - mean(y))
}
Z = center_apply(Z) #'recenter contols; this is mostly not needed
Delta=1
bracketed_wage<-bracket(mydata$lnw,Delta=Delta,1,5)
bracketed_wage$yM<-(bracketed_wage$yL+bracketed_wage$yU)/2
mydata$yL<-bracketed_wage$yL
mydata$yU<-bracketed_wage$yU
mydata$yM<-bracketed_wage$yM


### 


y_name   <- "yL";
d_name    <- "female";
form_z    <- "-1 + (exp1+exp2+exp3+exp4)*(educ+occ2+ind2+ms+reg)+(educ+occ2+ind2+ms+reg)"


a_p=first_stage_lasso(data=mydata,y_name="lnw",d_name="female",form_z=form_z)
##
b.hat.fit<-lm(y_residual~treat_residual-1,data=a_p)
b.hat<-coef(b.hat.fit)
se.hat<-sqrt(diag(vcovHC(b.hat.fit, type = 'HC')))
CI<-c(b.hat-1.96*se.hat, b.hat+1.96*se.hat)
res_table[1,1]<-b.hat
res_table[1,3:4]<-CI


a_M=first_stage_lasso(data=mydata,y_name="yM",d_name="female",form_z=form_z)
bM.hat.fit<-lm(y_residual~treat_residual-1,data=a_M)
bM.hat<-coef(bM.hat.fit)
seM.hat<-sqrt(diag(vcovHC(bM.hat.fit, type = 'HC')))
CI.M<-c(bM.hat-1.96*seM.hat, bM.hat+1.96*seM.hat)
res_table[2,1]<-bM.hat
res_table[2,3:4]<-CI.M



a_UL=first_stage_lasso(data=mydata,y_name="yL",d_name="female",form_z=form_z)
a_UL$yUBG<-a_UL$y_residual + Delta*(a_UL$treat_residual>0)-Delta/2 
bU.hat.fit<-lm(yUBG~treat_residual-1,data=a_UL)
bU.hat<-coef(bU.hat.fit)
seU.hat<-sqrt(diag(vcovHC(bU.hat.fit, type = 'HC')))





a_UL$yLBG<-a_UL$y_residual + Delta*(a_UL$treat_residual<0)-Delta/2 
bL.hat.fit<-lm(yLBG~treat_residual-1,data=a_UL)
bL.hat<-coef(bL.hat.fit)
seL.hat<-sqrt(diag(vcovHC(bL.hat.fit, type = 'HC')))

CR<-c(bL.hat-1.96*seL.hat, bU.hat+1.96*seU.hat)

res_table[3,1:2]<-c(bL.hat,bU.hat)
res_table[3,3:4]<-CR


### complex lasso specification
mydata$Sign<-a_UL$treat_residual>0
riesz.fit<-rlassologit(as.formula(paste0("Sign~",form_z )),data=mydata)
riesz_fittedvalue<-Delta*predict(riesz.fit,Z,type="response")
a_UL$yUBG<-a_UL$y_residual + Delta*(a_UL$treat_residual>0)-riesz_fittedvalue


bU.hat.fit<-lm(yUBG~treat_residual-1,data=a_UL)
bU.hat.complex<-coef(bU.hat.fit)
seU.hat.complex<-sqrt(diag(vcovHC(bU.hat.fit, type = 'HC')))

riesz_fittedvalue<-Delta*(1-predict(riesz.fit,Z,type="response"))
a_UL$yLBG<-a_UL$y_residual + Delta*(a_UL$treat_residual<0)-riesz_fittedvalue 
bL.hat.fit<-lm(yLBG~treat_residual-1,data=a_UL)
bL.hat.complex<-coef(bL.hat.fit)
seL.hat.complex<-sqrt(diag(vcovHC(bL.hat.fit, type = 'HC')))
CR.complex<-c(bL.hat.complex-1.96*seL.hat.complex, bU.hat.complex+1.96*seU.hat.complex)

res_table[4,1:2]<-c(bL.hat.complex,bU.hat.complex)
res_table[4,3:4]<-CR.complex
### save my results 

rownames(res_table)<-c("$Y$", "$Y_M$", "Simple","Complex")
tab<- xtable(res_table, digits=3)	
print.xtable(tab,type="latex", paste0(directoryname,"/results/Lasso1.txt"))