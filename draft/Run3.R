library(hdm)
library(randomForest)
library(ranger)
library(xtable)

directoryname<-"/Users/virasemenova/Dropbox (MIT)/JMP_Submission/JoE/code/"
setwd(directoryname)

source("empirical_application/Functions.R")
mydata<-read.csv("data/cps2015.csv")

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

### first stage ###

y_name   <- "yL";
d_name    <- "female";
form_z    <- "-1 + (exp1+exp2+exp3+exp4)*(educ+occ2+ind2+ms+reg)+(educ+occ2+ind2+ms+reg)"


a_p=first_stage_rf(Y=mydata$lnw,D=mydata$female,Z=Z)
##
b.hat.fit<-lm(y_residual~treat_residual-1,data=a_p)
b.hat<-coef(b.hat.fit)


a_M=first_stage_rf(Y=mydata$yM,D=mydata$female,Z=Z)
bM.hat.fit<-lm(y_residual~treat_residual-1,data=a_M)
bM.hat<-coef(bM.hat.fit)

a_UL=first_stage_rf(Y=mydata$yM,D=mydata$female,Z=Z)
a_UL$yUBG<-a_UL$y_residual + Delta*(a_UL$treat_residual>0)-Delta/2 
bU.hat.fit<-lm(yUBG~treat_residual-1,data=a_UL)
bU.hat<-coef(bU.hat.fit)

a_UL$yLBG<-a_UL$y_residual + Delta*(a_UL$treat_residual<0)-Delta/2 
bL.hat.fit<-lm(yLBG~treat_residual-1,data=a_UL)
bL.hat<-coef(bL.hat.fit)

##### random forest option


a_p=first_stage_lasso(y_name="lnw",data=mydata)
##
b.hat.fit<-lm(y_residual~treat_residual-1,data=a_p)
b.hat<-coef(b.hat.fit)


a_M=first_stage_lasso(y_name="yM",data=mydata)
bM.hat.fit<-lm(y_residual~treat_residual-1,data=a_M)
bM.hat<-coef(bM.hat.fit)

a_UL=first_stage_lasso(y_name="yL",data=mydata)
a_UL$yUBG<-a_UL$y_residual + Delta*(a_UL$treat_residual>0)-Delta/2 
bU.hat.fit<-lm(yUBG~treat_residual-1,data=a_UL)
bU.hat<-coef(bU.hat.fit)

a_UL$yLBG<-a_UL$y_residual + Delta*(a_UL$treat_residual<0)-Delta/2 
bL.hat.fit<-lm(yLBG~treat_residual-1,data=a_UL)
bL.hat<-coef(bL.hat.fit)