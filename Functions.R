### bracket the outcome
bracket<-function(y,Delta,minb=-8, maxb=8,...) {
  yL<-y
  yU<-y
  brackets<-seq(minb,maxb,by=Delta)
  sample_size<-length(y)
  for (i in 1:sample_size) {
    for (j in 1:(length(brackets)-1)) {
      if (brackets[j] <=y[i]  && y[i]< brackets[j+1]) {
        yL[i]<-brackets[j]
        yU[i]<-brackets[j+1]
      }
    }
  }
  return(data.frame(yL=yL,yU=yU,y=y))
}

### the function max(x,0)

vabs<-function(x){
  return(x*as.numeric(x>0))
}

## OLS

OLS<-function(covariate,outcome,weights=NULL) {
  outcome<-as.matrix(outcome)
  b.hat<-matrix(0,dim(covariate)[2],dim(outcome)[2])
  if (is.null(weights)) {
    weights<-rep(1,dim(covariate)[1])
  }
 
  for (j in 1:dim(outcome)[2]) {
    b.hat[,j]<-solve(cov.wt(covariate,wt=weights)$cov)%*%apply(covariate*matrix(rep(outcome[,j],dim(covariate)[2]),ncol=dim(covariate)[2]),2,weighted.mean,weights)
  }
  return(b.hat)
}

## simulate data


simulate_data_2d<-function(seed,decay,sample_size,p,beta,cons.x,cons.eta,Sigma,Delta,rho,control.name) {
  ### Generate data
  
  set.seed(seed)
  ### Generate controls
  gamma<-1/(1:p)^decay
  z<-matrix(rnorm(sample_size*p),ncol=p)
  firstrow<-rho^c(0:(p-1))
  Omega<-toeplitz(firstrow)
  z<- z%*%chol(Omega)
  V<-matrix(rnorm(2*sample_size),ncol=2)
  V<-V%*%chol(Sigma)
  
  x1<-(z%*%gamma)*cons.x[1]+V[,1]
  x2<-(z%*%gamma)*cons.x[2]+V[,2]
  
  
  x<-cbind(x1,x2)
  # Outcome
  eta<-1/(1:p)^decay
  #eta<-rep(0,p)
  xi<-rnorm(sample_size)
  y<- x%*%beta + (z%*%eta)*cons.eta + xi
  
  ### bracket the outcome
  datay<-bracket(y,Delta)
 
  
  data<-data.frame(X1=x1,X2=x2,V1=V[,1],V2=V[,2],y=datay$y,yL=datay$yL, yU=datay$yU,z=z)
  
  colnames(data)[8:length(colnames(data))]<-control.name
  return(data)
}



estimate_sigma<-function(data,num_circle,sample_size,Delta,treat.name,outcome.name,seed=NULL,fitted_outcome_name=NULL,partial_out=FALSE,weighted=FALSE, ...) {
  
  
  
  
  if (is.null(seed)) {
    inds<-1:sample_size
    weights<-rep(1,sample_size)
   
  } else {
    set.seed(seed)
    if (weighted==FALSE) {
      inds<-sample(1:sample_size,sample_size,replace=TRUE)
      weights<-rep(1,sample_size)
    } else {
      inds<-1:sample_size
      weights<-rexp(sample_size)
      weights<-weights/mean(weights)
    }
    
  }
  
  
  covariate<-data[inds,treat.name]
  outcomeL<-as.matrix(data[inds,outcome.name])
  outcomeU<-outcomeL+Delta
  
  ## support function estimate
  sigma.hat<-matrix(0,num_circle,dim(outcomeL)[2])
  
  ### covariance matrix
  Sigma.hat<-cov.wt(covariate,wt=weights)$cov
  for (k in 1:num_circle) {
    #print(k)
    alpha = 2*pi/num_circle*k
    q = c(cos(alpha),sin(alpha))
    
    #print(Sigma.hat)
    
    direction = t(q)%*%solve(Sigma.hat)%*%t(covariate)
    ### discrete version
    yubg = outcomeL
    yubg[direction>0,] =  outcomeU[direction>0,]
   
    if (partial_out==TRUE) {
      ## partial out the controls (see Remark 4.2)
      gammaL.X<-as.matrix(data[inds,fitted_outcome_name])
      outcome = yubg-gammaL.X-Delta/2 
    } else {
      outcome = yubg
    }
    if (weighted==TRUE) {
      b.hat=OLS(covariate,outcome,weights=weights) 
    } else {
      b.hat=OLS(covariate,outcome)
    }
    
    
    sigma.hat[k,] = t(q)%*% b.hat
  }
  return(sigma.hat)
}



#### Chandrasekhar, Chernozhukov, Molinari, Schrimpf (2012) et al only

first_stage_ols<-function(covariate, outcome,weights) {
  b.hat=OLS(covariate=as.matrix(covariate),outcome=as.matrix(outcome),weights=weights)
  fitted_value = as.matrix(covariate)%*%b.hat
  residual = outcome-fitted_value
  return(residual)
}

estimate_sigma_series<-function(data,num_circle,sample_size,Delta,treat.name,outcome.name,control.name,seed=NULL,fitted_outcome_name=NULL,partial_out=FALSE,weighted=FALSE, ...) {
  
  
  
  
  if (is.null(seed)) {
    inds<-1:sample_size
    weights<-rep(1,sample_size)
  } else {
    set.seed(seed)
    
    if (weighted == TRUE) {
      inds<-1:sample_size
      weights<-rexp(sample_size)
    } else {
      inds<-sample(1:sample_size,sample_size,replace=TRUE)
      weights<-rep(1,sample_size)
    }
  }
  
  
  covariate<-data[inds,treat.name]
  outcomeL<-as.matrix(data[inds,outcome.name])
  outcomeU<-outcomeL+Delta
  
  fitted_residuals<-first_stage_ols(data[inds,control.name],data[inds,treat.name],weights=weights)
  
  ## support function estimate
  sigma.hat<-matrix(0,num_circle,dim(outcomeL)[2])
  
  ### covariance matrix
  Sigma.hat<-cov.wt(covariate,wt=weights)$cov
  for (k in 1:num_circle) {
    #print(k)
    alpha = 2*pi/num_circle*k
    q = c(cos(alpha),sin(alpha))
    
    #print(Sigma.hat)
    
    direction = t(q)%*%solve(Sigma.hat)%*%t(fitted_residuals)
    ### discrete version
    yubg = outcomeL
    yubg[direction>0,] =  outcomeU[direction>0,]
    
    if (partial_out==TRUE) {
      ## partial out the controls (see Remark 4.2)
      
      cov<-as.matrix(data[inds,control.name])
      
      b.hat=OLS(covariate=cov,outcome=outcomeL,weights=weights)
      gammaL.X<-as.matrix(cov)%*%b.hat
      outcome = yubg-gammaL.X-Delta/2 
    } else {
      outcome = yubg
    }
    b.hat=OLS(fitted_residuals,outcome,weights=weights) 
    sigma.hat[k,] = t(q)%*% b.hat
  }
  return(sigma.hat)
}


### this function is not used anywhere in the 2d-design 
### coding upper bound estimate for 1d design

estimate_upper_bound<-function(data,sample_size,treat.name,outcome.name,fitted_outcome_name=NULL,fitted_delta_name=NULL,partial_out=FALSE) {
  
  if (is.null(seed)) {
    inds<-1:sample_size
    
  } else {
    set.seed(seed)
    inds<-sample(1:sample_size,sample_size,replace=TRUE)
  }
  
  ### retrieve data and fitted values
  covariate<-data[inds,treat.name]
  outcomeL<-as.matrix(data[inds,outcomeL.name])
  outcomeU<-as.matrix(data[inds,outcomeU.name])
  gammaL.X<-as.matrix(data[inds,fitted_outcome_name])
  Delta.X<-as.matrix(data[inds,fitted_delta_name])
  
  
  direction = covariate>0
  ### estimate best-case outcome
  yubg = outcomeL
  yubg[direction>0,] =  outcomeU[direction>0,]
  ## best-case outcome's residual
  outcome = yubg-gammaL.X-Delta.X/2 
  
  b.hat=OLS(covariate,outcome) 
  return(b.hat)
}
