##  This are the functions used in the Metropolis Hasting algorithm (MHA.R file)

######################################################
# Likelihood function
likelihood = function(Calibr_parameters){
  
  source("R calibration\\Main_model_calibration.R") # run the model
  
  output <- f.calibr.output(results_no_screen) # extract the outputs
  
  Predict <- cbind(output$rate_outcomes_m, output$rate_outcomes_f) #data prediction with given parameters
  
  likelihood <- (f.GOF.calc(1, m.GOF, Targets, SE, Predict))[1,]   # log likelihood instead of sum of squared errors
  
  sum_likelihoods = sum(likelihood)
  
  return(sum_likelihoods)

}

###########################################################
# Prior function


prior_unif = function(fitted_params){ #uniform distribution for priors
  
  Param_prior <- rep(0, nrow(fitted_params))
  
  for(i in 1:nrow(fitted_params)){ 
    if(fitted_params[i,1] >0){
      Param_prior[i] <- dunif(fitted_params[i,1], min=fitted_params[i,1]*0.9, max=fitted_params[i,1]*1.1, log = T)
    } else{
      Param_prior[i] <- dunif(fitted_params[i,1], min=fitted_params[i,1]*1.1, max=fitted_params[i,1]*0.9, log = T)
    }
  }
  return(sum(Param_prior))
  
}


#############################################################

posterior = function(param){
  return (likelihood(param) + prior_unif(param))
}


# The Distributions proposed:
# 4 = Uniform for parameters on age and sex onset and shape for time since onset
# 5 - Normal for C.age.80.undiag.mort
# 7 = Dirichlet for probabilities
# 8 = Truncated normal


####################################################################################
f.proposal.param = function(param){
  
  # Sample for parameters with the uniform distribution
  param[which(param[ ,"Distribution"]==4), "Param1"] <-  param[which(param[ ,"Distribution"]==4),"Mean"] - param[which(param[ ,"Distribution"]==4),"step"]*param[which(param[ ,"Distribution"]==4),"Mean"]
  param[which(param[ ,"Distribution"]==4), "Param2"] <-  param[which(param[ ,"Distribution"]==4),"Mean"] + param[which(param[ ,"Distribution"]==4),"step"]*param[which(param[ ,"Distribution"]==4),"Mean"]
  
  
  #return(rnorm(length(param), mean = param, sd= step))
  
  
  # Re-sample parameters if the distribution is out of the plausible values
  param["P.onset_age", "Param1"] <- replace(param["P.onset_age", "Param1"], param["P.onset_age", "Param1"] <1, 1)
  param["RR.onset_sex", "Param1"] <- replace(param["RR.onset_sex", "Param1"], param["RR.onset_sex", "Param1"] <1, 1)
  param["P.sympt.diag_LGBC", "Param1"] <- replace(param["P.sympt.diag_LGBC", "Param1"], param["P.sympt.diag_LGBC", "Param1"] <5.000000e-02, 5.000000e-02)
  param["P.sympt.diag_A_HGBC", "Param1"] <- replace(param["P.sympt.diag_A_HGBC", "Param1"], param["P.sympt.diag_A_HGBC", "Param1"] < 1.000000e-02, 1.000000e-02)
  param["P.sympt.diag_B_HGBC", "Param1"] <- replace(param["P.sympt.diag_B_HGBC", "Param1"], param["P.sympt.diag_B_HGBC", "Param1"] <1.010000e+00, 1.010000e+00)
  param["P.sympt.diag_Age80_HGBC", "Param1"] <- replace(param["P.sympt.diag_Age80_HGBC", "Param1"], param["P.sympt.diag_Age80_HGBC", "Param1"] < 5.000000e-01, 5.000000e-01)
  param["shape.t.StI.StII", "Param1"] <- replace(param["shape.t.StI.StII", "Param1"], param["shape.t.StI.StII", "Param1"] < 1.001000e+00, 1.001000e+00)
  param["shape.t.StII.StIII", "Param1"] <- replace(param["shape.t.StII.StIII", "Param1"], param["shape.t.StII.StIII", "Param1"] < 1.001000e+00, 1.001000e+00)
  param["shape.t.StIII.StIV", "Param1"] <- replace(param["shape.t.StIII.StIV", "Param1"], param["shape.t.StIII.StIV", "Param1"] < 1.001000e+00, 1.001000e+00)
  param["P.LGtoHGBC", "Param1"] <- replace(param["P.LGtoHGBC", "Param1"], param["P.LGtoHGBC", "Param1"] < 0.001, 0.001)
  
  param["P.onset_age", "Param2"] <- replace(param["P.onset_age", "Param2"], param["P.onset_age", "Param2"] > 1.10, 1.10)
  param["RR.onset_sex", "Param2"] <- replace(param["RR.onset_sex", "Param2"], param["RR.onset_sex", "Param2"] > 1.5, 1.5)
  param["P.sympt.diag_LGBC", "Param2"] <- replace(param["P.sympt.diag_LGBC", "Param2"], param["P.sympt.diag_LGBC", "Param2"] >0.10, 0.10)
  param["P.sympt.diag_A_HGBC", "Param2"] <- replace(param["P.sympt.diag_A_HGBC", "Param2"], param["P.sympt.diag_A_HGBC", "Param2"] >0.20, 0.20)
  param["P.sympt.diag_B_HGBC", "Param2"] <- replace(param["P.sympt.diag_B_HGBC", "Param2"], param["P.sympt.diag_B_HGBC", "Param2"] > 1.70, 1.70)
  param["P.sympt.diag_Age80_HGBC", "Param2"] <- replace(param["P.sympt.diag_Age80_HGBC", "Param2"], param["P.sympt.diag_Age80_HGBC", "Param2"] >1.00000e+00, 1.00000e+00)
  param["shape.t.StI.StII", "Param2"] <- replace(param["shape.t.StI.StII", "Param2"], param["shape.t.StI.StII", "Param2"] > 1.2, 1.2)
  param["shape.t.StII.StIII", "Param2"] <- replace(param["shape.t.StII.StIII", "Param2"], param["shape.t.StII.StIII", "Param2"] > 1.2, 1.2)
  param["shape.t.StIII.StIV", "Param2"] <- replace(param["shape.t.StIII.StIV", "Param2"], param["shape.t.StIII.StIV", "Param2"] > 1.2, 1.2)
  param["P.LGtoHGBC", "Param2"] <- replace(param["P.LGtoHGBC", "Param2"], param["P.LGtoHGBC", "Param2"] >0.1, 0.1)
  

  # Re-calculate Param1 and Param2 based on the means and the distributions for beta distribution
  
  # Update alpha and beta for parameters with beta distribution
  #param["P.onset", "Mean"] <- truncnorm::rtruncnorm(1, a=0, b=1, mean=param["P.onset", "Mean"], sd=param["P.onset","step"])
  #param["P.onset_low.risk", "Mean"] <- truncnorm::rtruncnorm(1, a=0, b=1, mean=param["P.onset_low.risk", "Mean"], sd=param["P.onset_low.risk","step"])
  #param["P.LGtoHGBC", "Mean"] <- truncnorm::rtruncnorm(1, a=0.001, b=0.1, mean=param["P.LGtoHGBC", "Mean"], sd=param["P.LGtoHGBC","step"])
  
  #param["P.onset", "Param1"]<- ((1-param["P.onset", "Mean"])/(0.01^2)^2 -1/param["P.onset", "Mean"])*param["P.onset", "Mean"]^2 #Use SD reported in the studies for prior calc
  #param["P.onset", "Param2"]<- param["P.onset", "Param1"]*(1/param["P.onset", "Mean"]-1)
  
  #param["P.onset_low.risk", "Param1"]<- ((1-param["P.onset_low.risk", "Mean"])/(0.05^2)^2 -1/param["P.onset_low.risk", "Mean"])*param["P.onset_low.risk", "Mean"]^2
  #param["P.onset_low.risk", "Param2"]<- param["P.onset_low.risk", "Param1"]*(1/param["P.onset_low.risk", "Mean"]-1)
  
  #param["P.LGtoHGBC", "Param1"]<- ((1-param["P.LGtoHGBC", "Mean"])/(0.05^2)^2 -1/param["P.LGtoHGBC", "Mean"])*param["P.LGtoHGBC", "Mean"]^2 #Use 10% SD 
 # param["P.LGtoHGBC", "Param2"]<- param["P.LGtoHGBC", "Param1"]*(1/param["P.LGtoHGBC", "Mean"]-1)
  
   #proposal <- runif(1, min = param - step*param, max = param + step*param)
  
  
  return(param)
  
}


##################################################################################################################

# DRAFT FUNCTIONS


###########################################################


proposalfunction = function(param, Param_sets, N_sets){
  
  # Sample for parameters with the uniform distribution
  Param_sets[ ,"shape.t.StI.StII"] <<-  runif(N_sets, min = param["shape.t.StI.StII","Mean"] - param["shape.t.StI.StII","step"]*param["shape.t.StI.StII","Mean"], max = param["shape.t.StI.StII","Mean"] + param["shape.t.StI.StII","step"]*param["shape.t.StI.StII","Mean"])
  Param_sets[ ,"shape.t.StII.StIII"] <<-  runif(N_sets, min = param["shape.t.StII.StIII","Mean"] - param["shape.t.StII.StIII","step"]*param["shape.t.StII.StIII","Mean"], max = param["shape.t.StII.StIII","Mean"] + param["shape.t.StII.StIII","step"]*param["shape.t.StII.StIII","Mean"])
  Param_sets[ ,"shape.t.StIII.StIV"] <<-  runif(N_sets, min = param["shape.t.StIII.StIV","Mean"] - param["shape.t.StIII.StIV","step"]*param["shape.t.StIII.StIV","Mean"], max = param["shape.t.StIII.StIV","Mean"] + param["shape.t.StIII.StIV","step"]*param["shape.t.StIII.StIV","Mean"])
  Param_sets[ ,"C.age.80.undiag.mort"] <<-  runif(N_sets, min = param["C.age.80.undiag.mort","Mean"] + param["C.age.80.undiag.mort","step"]*param["C.age.80.undiag.mort","Mean"], max = param["C.age.80.undiag.mort","Mean"] - param["C.age.80.undiag.mort","step"]*param["C.age.80.undiag.mort","Mean"])
  
  # Sample for parameters with the truncated normal distribution
  Param_sets[ ,"P.onset_age"] <<-  truncnorm::rtruncnorm(N_sets, a =Calibr_parameters["P.onset_age", "Param1"], b=Calibr_parameters["P.onset_age", "Param2"], 
                                                         mean=Calibr_parameters["P.onset_age", "Mean"], sd=Calibr_parameters["P.onset_age","step"]) #re-set the sd 
  
  Param_sets[ ,"RR.onset_sex"] <<-  truncnorm::rtruncnorm(N_sets, a =Calibr_parameters["RR.onset_sex", "Param1"], b=Calibr_parameters["RR.onset_sex", "Param2"], 
                                                          mean=Calibr_parameters["RR.onset_sex", "Mean"], sd=Calibr_parameters["RR.onset_sex","step"]) #re-set the sd 
  
  Param_sets[ ,"P.sympt.diag_LGBC"] <<-  truncnorm::rtruncnorm(N_sets, a =Calibr_parameters["P.sympt.diag_LGBC", "Param1"], b=Calibr_parameters["P.sympt.diag_LGBC", "Param2"], 
                                                               mean=Calibr_parameters["P.sympt.diag_LGBC", "Mean"], sd=Calibr_parameters["P.sympt.diag_LGBC","step"]) #re-set the sd 
  Param_sets[ ,"P.sympt.diag_A_HGBC"] <<-  truncnorm::rtruncnorm(N_sets, a =Calibr_parameters["P.sympt.diag_A_HGBC", "Param1"], b=Calibr_parameters["P.sympt.diag_A_HGBC", "Param2"], 
                                                                 mean=Calibr_parameters["P.sympt.diag_A_HGBC", "Mean"], sd=Calibr_parameters["P.sympt.diag_A_HGBC","step"]) #re-set the sd 
  Param_sets[ ,"P.sympt.diag_Age80_HGBC"] <<-  truncnorm::rtruncnorm(N_sets, a =Calibr_parameters["P.sympt.diag_Age80_HGBC", "Param1"], b=Calibr_parameters["P.sympt.diag_Age80_HGBC", "Param2"], 
                                                                     mean=Calibr_parameters["P.sympt.diag_Age80_HGBC", "Mean"], sd=Calibr_parameters["P.sympt.diag_Age80_HGBC","step"]) #re-set the sd 
  
  
  # Re-calculate Param1 and Param2 based on the means and the distributions
  
  # Update alpha and beta for parameters with beta distribution
  Calibr_parameters["P.onset", "Mean"] <- truncnorm::rtruncnorm(1, a=0, b=1, mean=param["P.onset", "Mean"], sd=param["P.onset","step"])
  Calibr_parameters["P.onset_low.risk", "Mean"] <- truncnorm::rtruncnorm(1, a=0, b=1, mean=param["P.onset_low.risk", "Mean"], sd=param["P.onset_low.risk","step"])
  
  Calibr_parameters["P.onset", "Param1"]<- ((1-Calibr_parameters["P.onset", "Mean"])/(0.01^2)^2 -1/Calibr_parameters["P.onset", "Mean"])*Calibr_parameters["P.onset", "Mean"]^2
  Calibr_parameters["P.onset", "Param2"]<- Calibr_parameters["P.onset", "Param1"]*(1/Calibr_parameters["P.onset", "Mean"]-1)
  
  Calibr_parameters["P.onset_low.risk", "Param1"]<- ((1-Calibr_parameters["P.onset_low.risk", "Mean"])/(0.01^2)^2 -1/Calibr_parameters["P.onset_low.risk", "Mean"])*Calibr_parameters["P.onset_low.risk", "Mean"]^2
  Calibr_parameters["P.onset_low.risk", "Param2"]<- Calibr_parameters["P.onset_low.risk", "Param1"]*(1/Calibr_parameters["P.onset_low.risk", "Mean"]-1)
  
  Param_sets[ ,"P.onset"] <<- rbeta(N_sets, shape1=param["P.onset", "Param1"], shape2=param["P.onset", "Param2"])
  Param_sets[ ,"P.onset_low.risk"] <<- rbeta(N_sets, shape1=param["P.onset_low.risk", "Param1"], shape2=param["P.onset_low.risk", "Param2"])
  
  #proposal <- runif(1, min = param - step*param, max = param + step*param)
  
  #return(rnorm(length(param), mean = param, sd= step))
  
  return(Param_sets)
  
}


###############################################
# Analyse the output of the MHA
f.analyse.param.MHA <- function(){
  
  P.onset <- cbind(Param1[ , "P.onset"], Param2[ , "P.onset"])
  P.onset_low.risk <- cbind(Param1[ , "P.onset_low.risk"], Param2[ , "P.onset_low.risk"])
  P.onset_age <- cbind(Param1[ , "P.onset_age"], Param2[ , "P.onset_age"])
  RR.onset_sex <- cbind(Param1[ , "RR.onset_sex"], Param2[ , "RR.onset_sex"])
  
  P.sympt.diag_LGBC <- cbind(Param1[ , "P.sympt.diag_LGBC"], Param2[ , "P.sympt.diag_LGBC"])
  P.sympt.diag_A_HGBC <- cbind(Param1[ , "P.sympt.diag_A_HGBC"], Param2[ , "P.sympt.diag_A_HGBC"])
  P.sympt.diag_B_HGBC <- cbind(Param1[ , "P.sympt.diag_B_HGBC"], Param2[ , "P.sympt.diag_B_HGBC"])
  P.sympt.diag_Age80_HGBC <- cbind(Param1[ , "P.sympt.diag_Age80_HGBC"], Param2[ , "P.sympt.diag_Age80_HGBC"])
  shape.t.StI.StII <- cbind(Param1[ , "shape.t.StI.StII"], Param2[ , "shape.t.StI.StII"])
  shape.t.StII.StIII <- cbind(Param1[ , "shape.t.StII.StIII"], Param2[ , "shape.t.StII.StIII"])
  shape.t.StIII.StIV <- cbind(Param1[ , "shape.t.StIII.StIV"], Param2[ , "shape.t.StIII.StIV"])
  
  output_params <- list(P.onset, P.onset_low.risk, P.onset_age, RR.onset_sex,
                        P.sympt.diag_LGBC, P.sympt.diag_A_HGBC, P.sympt.diag_B_HGBC,
                        P.sympt.diag_Age80_HGBC, shape.t.StI.StII, shape.t.StII.StIII,
                        shape.t.StIII.StIV)
}
