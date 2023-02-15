# drafts for Bayesian claibration functions

#####################################################################

f.generate_parameters_Bayes <- function(Param_sets, Calibr_parameters, N_sets){
  
  # 1 - beta, 4 -uniform, 5 - normal, 7 -truncated normal
  
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"P.onset")
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"P.onset_low.risk")
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"P.onset_age")
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"RR.onset_sex")
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"P.sympt.diag_LGBC")
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"P.sympt.diag_B_HGBC")
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"P.sympt.diag_Age80_HGBC")
  #f.update.dist(N_sets, Calibr_parameters, Param_sets,"C.age.80.undiag.mort")
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"shape.t.StI.StII")
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"shape.t.StII.StIII")
  f.update.dist(N_sets, Calibr_parameters, Param_sets,"shape.t.StIII.StIV")
  
  Param_sets[ ,"C.age.80.undiag.mort"] <-  runif(N_sets, min=Calibr_parameters["C.age.80.undiag.mort", "Param2"], max=Calibr_parameters["C.age.80.undiag.mort", "Param1"])
  
  Param_sets
}


f.generate_parameters_Bayes(Param_sets, Calibr_parameters, N_sets)

############################################################################  
## Function to update the distributions in the MHA

f.update.dist <- function(N_sets, Calibr_parameters, replaced_params, name_param){
  
  if(Calibr_parameters[name_param, "Distribution"]==1){
    replaced_params[1:N_sets,name_param] <- rbeta(N_sets, shape1=Calibr_parameters[name_param, "Param1"], shape2=Calibr_parameters[name_param, "Param2"])
    
  } else if(Calibr_parameters[name_param, "Distribution"]==4){
    replaced_params[1:N_sets,name_param] <- runif(N_sets, min=Calibr_parameters[name_param, "Param1"]- Calibr_parameters[name_param,"step"]*Calibr_parameters[ ,"Mean"], max=Calibr_parameters[name_param, "Param1"] + Calibr_parameters[name_param,"step"]*Calibr_parameters[ ,"Mean"])
    
  } else if(Calibr_parameters[name_param, "Distribution"]==5){
    replaced_params[1:N_sets,name_param] <- rnorm(N_sets, mean = Calibr_parameters[name_param, "Param1"], sd = Calibr_parameters[name_param,"step"])
    
  } else if(Calibr_parameters[name_param, "Distribution"]==7){
    replaced_params[1:N_sets,name_param] <- truncnorm::rtruncnorm(N_sets, a =Calibr_parameters[name_param, "Param1"], b=Calibr_parameters[name_param, "Param2"], 
                                                                  mean=Calibr_parameters[name_param, "Mean"], sd=Calibr_parameters[name_param,"step"]) #re-set the sd 
  }
  return(replaced_params)
}