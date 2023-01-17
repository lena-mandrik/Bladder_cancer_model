##  This is the code to run calibration either with a Metropolis Hasting algorithm

# Use randomly calibrated parameters as priors
fitted_params <- (read.table("R calibration/Outputs/parameters_fit.txt", header = F, row.names=1)) #exclude the parameter on RR death for no smokers, as it was calculated
Calibr_parameters[,1] <- fitted_params[,1]

sd_params <- abs(fitted_params[,1]*0.1)

step <- abs(fitted_params[,1]*0.1) #10% step

######################################################
# Likelihood function
likelihood = function(Calibr_parameters){
  
  #Calibr_parameters[,1] <- param
  
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



###########################################################
proposalfunction = function(param){
  return(rnorm(length(param), mean = param, sd= step))
}


######## Metropolis algorithm ################


run_metropolis_MCMC = function(startvalue, iterations, Calibr_parameters){
  
  chain = array(dim = c(iterations+1,nrow(startvalue)))
  
  report_prob = matrix(nrow=iterations+1, ncol=1)
  
  posterior_rep =matrix(nrow=iterations+1, ncol=2)
  
  chain[1,] = t(startvalue[,1])
  
  for (i in 1:iterations){
    
    proposal = proposalfunction(t(chain[i,]))
    
    Calibr_parameters[,1] <<- proposal #save to the Calibr_parameters in Global env
    posterior_proposal = posterior(Calibr_parameters)
    
    Calibr_parameters[,1] <<- chain[i,] #save to the Calibr_parameters in Global env
    posterior_current = posterior(Calibr_parameters)
    
    probab = exp(posterior_proposal - posterior_current)
    
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    
    report_prob[i]=probab
    posterior_rep[i,1]=posterior_proposal
    posterior_rep[i,2]=posterior_current
    
  }
  
 # write.csv(report_prob, file="R calibration\\Outputs\\test_prob.csv")

  output <- cbind(chain,report_prob,posterior_rep)
  
  return(output)
}

####################################################################################################################


test <- run_metropolis_MCMC(fitted_params, 10, Calibr_parameters)











Distributions <-matrix(NA, nrow=nrow(fitted_params), ncol =4)

Distributions[ ,1] <- fitted_params[,1]
rownames(Distributions) <- rownames(fitted_params)

# The Distributions proposed:
# 4 = Uniform for parameters on age and sex onset and shape for time since onset
# 5 - Normal for C.age.80.undiag.mort
# 7 = Dirichlet for probabilities 

Distributions[ ,2] <- rep(4, nrow(fitted_params))  # the distribution to consider for priors. Use uniform in the beginning. c(7,7,4,4,rep(7,4),5,7,rep(4,3)) #Define the distribution for the priors
Distributions[ ,3] <- Distributions[ ,1]*0.9 # set the uniform distribution 90% from the parameter value
Distributions[ ,4] <- Distributions[ ,1]*1.1 # set the uniform distribution 90% from the parameter value








Param_sets[which(Params[, "Distribution"] ==4),2:N_sets] <- runif(length(Param_sets[which(Params[, "Distribution"] ==4),1])*(N_sets-1), 
                                                                  min=Params[which(Params[, "Distribution"] ==4), "Param1"], 
                                                                  max=Params[which(Params[, "Distribution"] ==4), "Param2"])




 

# Example: plot the likelihood profile of the slope a
slopevalues = function(x){
  
  Calibr_parameters[,1] <- x
  
  return(likelihood(Calibr_parameters))}

slopelikelihoods = lapply(x, slopevalues )
plot (seq(3, 7, by=.05), slopelikelihoods , type="l", xlab = "values of slope parameter a", ylab = "Log likelihood")




# Set prior distributions
# Check fitted parameters from the previous calibration









# Set the parameters for claibration

  Distributions["P.onset", c(3:4)] <- 
  Distributions["P.onset_low.risk", c(3:4)] <-
  Distributions["P.onset_age", c(3:4)] <-
  Distributions["P.onset_sex", c(3:4)] <-
  Distributions["P.sympt.diag_LGBC", c(3:4)] <-
  Distributions["P.sympt.diag_A_HGBC", c(3:4)] <-
  Distributions["P.sympt.diag_B_HGBC", c(3:4)] <-
  Distributions["P.sympt.diag_Age80_HGBC", c(3:4)] <-
  Distributions["C.age.80.undiag.mort", c(3:4)] <-
  Distributions["RR.All.Death.no_smoke", c(3:4)] <-
  Distributions["shape.t.StI.StII", c(3:4)] <-
  Distributions["shape.t.StII.StIII", c(3:4)] <-
  Distributions["shape.t.StIII.StIV", c(3:4)] <-

  






startvalue = c(4,0,10)
chain = run_metropolis_MCMC(startvalue, 10000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
