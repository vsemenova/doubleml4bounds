first_stage<-function(data,
                      treat.name, outcome.name,control.name,
                      method.treat = gamlr,
                      method.outcome = gamlr,num_splits,INDS,...) {
  
  # create sample splitting pattern
  
  
  true_treat=data[,c("V1","V2")]
  true_outcome=data[,outcome.name]
  ### Residualize sales
  ### ... are settings of method.controls
  # show("Partialling out treatments ...")
  
  treat<-data[,treat.name]
  outcome<-data[,outcome.name]
  controls<-data[,control.name]
  
  t<-remove_wrapper(controls = controls,
                    target =treat,
                    method = method.treat,
                    num_splits=num_splits,INDS=INDS,...)
  treat<-t$res
  treat.error<-t$mse
  ### Residualize prices
  # show("Partialling out outcome ...")
  out<-remove_wrapper(controls = controls,
                      target = outcome,
                      method = method.outcome,
                      num_splits=num_splits,INDS=INDS,...)
  outcome<-out$res
  outcome.error<-out$mse
  
  fitted_outcome<-data[,outcome.name]-out$res
  
  return(data.frame(outcome = outcome, 
                    fitted_outcome=fitted_outcome,
                    true_outcome=true_outcome,
                    treat = treat,
                    true_treat=true_treat))
}