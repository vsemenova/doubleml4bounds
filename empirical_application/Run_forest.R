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
Delta=2
bracketed_wage<-bracket(mydata$lnw,Delta=Delta,1,5)
bracketed_wage$yM<-(bracketed_wage$yL+bracketed_wage$yU)/2
mydata$yL<-bracketed_wage$yL
mydata$yU<-bracketed_wage$yU
mydata$yM<-bracketed_wage$yM


### treatment residual
seed=1
D=mydata$female
Y=mydata$lnw
YM=mydata$yM
YL=mydata$yL

N<-length(D)
# Propensity score
s.hat<-rep(1,N)
y.hat<-rep(1,N)
yM.hat<-rep(1,N)
yL.hat<-rep(1,N)
## define sample splitting
set.seed(seed)
inds.train=sample(1:N,floor(N/2))
inds.eval=setdiff(1:N,inds.train)
form_z    <- "-1 + (exp1+exp2+exp3+exp4)+(educ+occ2+ind2+ms+reg)+(educ+occ2+ind2+ms+reg)"
print ("Estimate treatment probability, first half")
mydata$female<-as.factor(as.character(mydata$female))

fitted.rf.pscore<-ranger(as.formula(paste0("female~",form_z)),data=mydata[inds.train,],seed=1000,
                         probability=TRUE,respect.unordered.factors=TRUE)

s.hat[inds.eval]<-predict(fitted.rf.pscore,mydata[inds.eval,])$predictions[,2]

print ("Estimate treatment probability, second half")
fitted.rf.pscore<-ranger(as.formula(paste0("female~",form_z)),data=mydata[inds.eval,],seed=1000,
                         probability=TRUE,respect.unordered.factors=TRUE)
s.hat[inds.train]<-predict(fitted.rf.pscore,mydata[inds.train,])$predictions[,2]

fitted.rf.y<-ranger(as.formula(paste0("lnw~",form_z)),data=mydata[inds.train,],seed=1000,respect.unordered.factors=TRUE)
y.hat[inds.eval]<-predict(fitted.rf.y,mydata[inds.eval,])$predictions
fitted.rf.y<-ranger(as.formula(paste0("lnw~",form_z)),data=mydata[inds.eval,],seed=1000,respect.unordered.factors=TRUE)
y.hat[inds.train]<-predict(fitted.rf.y,mydata[inds.train,])$predictions

treat_residual<-as.numeric(as.character(mydata$female))-s.hat
y_residual<-mydata$lnw-y.hat

b.hat.fit<-lm(y_residual~treat_residual-1)
b.hat<-coef(b.hat.fit)
se.hat<-sqrt(diag(vcovHC(b.hat.fit, type = 'HC')))
CI<-c(b.hat-1.96*se.hat, b.hat+1.96*se.hat)
res_table[1,1]<-b.hat
res_table[1,3:4]<-CI

fitted.rf.y<-ranger(as.formula(paste0("yM~",form_z)),data=mydata[inds.train,],seed=1000,respect.unordered.factors=TRUE)
yM.hat[inds.eval]<-predict(fitted.rf.y,mydata[inds.eval,])$predictions
fitted.rf.y<-ranger(as.formula(paste0("yM~",form_z)),data=mydata[inds.eval,],seed=1000,respect.unordered.factors=TRUE)
yM.hat[inds.train]<-predict(fitted.rf.y,mydata[inds.train,])$predictions
yM_residual<-mydata$yM-yM.hat


bM.hat.fit<-lm(yM_residual~treat_residual-1)
bM.hat<-coef(bM.hat.fit)
seM.hat<-sqrt(diag(vcovHC(bM.hat.fit, type = 'HC')))
CI.M<-c(bM.hat-1.96*seM.hat, bM.hat+1.96*seM.hat)
res_table[2,1]<-bM.hat
res_table[2,3:4]<-CI.M

fitted.rf.y<-ranger(as.formula(paste0("yL~",form_z)),data=mydata[inds.train,],seed=1000,respect.unordered.factors=TRUE)
yL.hat[inds.eval]<-predict(fitted.rf.y,mydata[inds.eval,])$predictions
fitted.rf.y<-ranger(as.formula(paste0("yL~",form_z)),data=mydata[inds.eval,],seed=1000,respect.unordered.factors=TRUE)
yL.hat[inds.train]<-predict(fitted.rf.y,mydata[inds.train,])$predictions
yL_residual<-mydata$yL-yL.hat



yUBG<-yL_residual + Delta*(treat_residual>0)-Delta/2 
bU.hat.fit<-lm(yUBG~treat_residual-1)
bU.hat<-coef(bU.hat.fit)
seU.hat<-sqrt(diag(vcovHC(bU.hat.fit, type = 'HC')))





yLBG<-yL_residual + Delta*(treat_residual<0)-Delta/2 
bL.hat.fit<-lm(yLBG~treat_residual-1)
bL.hat<-coef(bL.hat.fit)
seL.hat<-sqrt(diag(vcovHC(bL.hat.fit, type = 'HC')))

CR<-c(bL.hat-1.96*seL.hat, bU.hat+1.96*seU.hat)

res_table[3,1:2]<-c(bL.hat,bU.hat)
res_table[3,3:4]<-CR





Sign<-treat_residual>0
riesz.fit<-rlassologit(as.formula(paste0("Sign~",form_z )),data=mydata)
riesz_fittedvalue<-Delta*predict(riesz.fit,Z,type="response")
yUBG<-yL_residual + Delta*(treat_residual>0)-riesz_fittedvalue


bU.hat.fit<-lm(yUBG~treat_residual-1)
bU.hat.complex<-coef(bU.hat.fit)
seU.hat.complex<-sqrt(diag(vcovHC(bU.hat.fit, type = 'HC')))

riesz_fittedvalue<-Delta*(1-predict(riesz.fit,Z,type="response"))
yLBG<-yL_residual + Delta*(treat_residual<0)-riesz_fittedvalue 
bL.hat.fit<-lm(yLBG~treat_residual-1)
bL.hat.complex<-coef(bL.hat.fit)
seL.hat.complex<-sqrt(diag(vcovHC(bL.hat.fit, type = 'HC')))
CR.complex<-c(bL.hat.complex-1.96*seL.hat.complex, bU.hat.complex+1.96*seU.hat.complex)

res_table[4,1:2]<-c(bL.hat.complex,bU.hat.complex)
res_table[4,3:4]<-CR.complex
### save my results 

rownames(res_table)<-c("$Y$", "$Y_M$", "Simple","Complex")
res_table<-apply(res_table,2,round,digits=3)
write.csv(res_table,paste0(directoryname,"/empirical_application/Tables/Forest.csv"))