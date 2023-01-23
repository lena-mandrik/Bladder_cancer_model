##  This is the code to run calibration either with a Metropolis Hasting algorithm

# Load starting parameters and set up the step
Calibr_parameters <- f.load.start.param(Calibr_parameters)

sd_params <- abs(fitted_params[,1]*0.1)


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


prior = function(param){ #uniform distribution for priors
  
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
  return (likelihood(param) + prior(param))
}


# The Distributions proposed:
# 4 = Uniform for parameters on age and sex onset and shape for time since onset
# 5 - Normal for C.age.80.undiag.mort
# 7 = Dirichlet for probabilities
# 8 = Truncated normal
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


f.proposal.param = function(param, N_sets){
  
  # Sample for parameters with the uniform distribution
  param[which(param[ ,"Distribution"]==4), "Param1"] <-  param[which(param[ ,"Distribution"]==4),"Mean"] - param[which(param[ ,"Distribution"]==4),"step"]*param[which(param[ ,"Distribution"]==4),"Mean"]
  param[which(param[ ,"Distribution"]==4), "Param2"] <-  param[which(param[ ,"Distribution"]==4),"Mean"] + param[which(param[ ,"Distribution"]==4),"step"]*param[which(param[ ,"Distribution"]==4),"Mean"]
  
    
  # Sample for parameters with the truncated normal distribution
  param[which(param[ ,"Distribution"]==7), "Param1"] <-   truncnorm::rtruncnorm(nrow(param[which(param[ ,"Distribution"]==7), ]), a=1, b=1.7, mean =param[which(param[ ,"Distribution"]==7), "Param1"], 
                                                               sd= param[which(param[ ,"Distribution"]==7), "Param2"])

  # Re-calculate Param1 and Param2 based on the means and the distributions for beta distribution
  
  # Update alpha and beta for parameters with beta distribution
  param["P.onset", "Mean"] <- truncnorm::rtruncnorm(1, a=0, b=1, mean=param["P.onset", "Mean"], sd=param["P.onset","step"])
  param["P.onset_low.risk", "Mean"] <- truncnorm::rtruncnorm(1, a=0, b=1, mean=param["P.onset_low.risk", "Mean"], sd=param["P.onset_low.risk","step"])
  
  param["P.onset", "Param1"]<- ((1-param["P.onset", "Mean"])/(0.01^2)^2 -1/param["P.onset", "Mean"])*param["P.onset", "Mean"]^2
  param["P.onset", "Param2"]<- param["P.onset", "Param1"]*(1/param["P.onset", "Mean"]-1)
  
  param["P.onset_low.risk", "Param1"]<- ((1-param["P.onset_low.risk", "Mean"])/(0.01^2)^2 -1/param["P.onset_low.risk", "Mean"])*param["P.onset_low.risk", "Mean"]^2
  param["P.onset_low.risk", "Param2"]<- param["P.onset_low.risk", "Param1"]*(1/param["P.onset_low.risk", "Mean"]-1)
  
   #proposal <- runif(1, min = param - step*param, max = param + step*param)
  
  #return(rnorm(length(param), mean = param, sd= step))
  
  return(param)
  
}


######## Metropolis algorithm ################

run_metropolis_MCMC = function(n_samples, Calibr_parameters, step){
  
  chain = array(dim = c(n_samples+1,nrow(Calibr_parameters)*2))
  
  report_prob = matrix(nrow=n_samples+1, 0, ncol=1)
  
  posterior_rep =matrix(nrow=n_samples+1, ncol=2)
  
  chain[1,] = t(c(Calibr_parameters[,3], Calibr_parameters[,4]))
  
  for (i in 1:iterations){
    
    #proposal = proposalfunction(t(chain[i,]), step)
    
    Calibr_parameters = f.proposal.param(Calibr_parameters, N_sets)
    chain[i,] <- t(c(Calibr_parameters[,"Param1"], Calibr_parameters[,"Param2"])) #record the sampled parameters
    
    Calibr_parameters[,1] <<- proposal #save to the Calibr_parameters in Global env
    posterior_proposal = posterior(Calibr_parameters)
    
    Calibr_parameters[,1] <<- chain[i,] #save to the Calibr_parameters in Global env
    posterior_current = posterior(Calibr_parameters)
    
    probab = exp(posterior_proposal - posterior_current)
    
    if (runif(1) < probab | runif(1) <0.05){ #accept all parameters with better fit and 5% with worse fit
      chain[i+1,] = proposal
      report_prob[i]=1
      
    }else{
      chain[i+1,] = chain[i,]
    }
    
    posterior_rep[i,1]=posterior_proposal
    posterior_rep[i,2]=posterior_current
    
  }
  
  if(i%%50==0){
    if(sum(report_prob[(i-100):i])/100 < 0.1){step = step*0.9}
  }
 # write.csv(report_prob, file="R calibration\\Outputs\\test_prob.csv")

  output <- cbind(chain,report_prob,posterior_rep)
  
  colnames(output) <- c(rownames(startvalue), "Accept", "posterior_proposal", "posterior_current")
  
  return(output)
}

####################################################################################################################


test <- run_metropolis_MCMC(fitted_params, 200, Calibr_parameters, step)


















