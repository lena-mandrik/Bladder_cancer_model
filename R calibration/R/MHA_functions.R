##  This are the functions used in the Metropolis Hasting algorithm (MHA.R file)

######################################################
# Likelihood function
likelihood = function(Calibr_parameters){
  
  source("R calibration\\Main_model_calibration.R") # run the model
  
  output <- f.calibr.output(results_no_screen) # extract the outputs
  
  Predict <- cbind(output$rate_outcomes_m, output$rate_outcomes_f) #data prediction with given parameters
  if(target_stat=="counts"){Predict <- as.matrix(sum_age_predictions(Predict, England.pop))}
 
  
  #likelihood <- f.GOF.MHA.calc(Targets, SE, Predict)  # log likelihood instead of sum of squared errors
  
  #sum_likelihoods = sum(likelihood)
  
  # Calculate the likelihood as the sum of log-likelihoods for each target variable
  log_likelihood <- 0
  
  for (i in 1:ncol(Targets)) {
    obs <- Targets[, i]
    sd <- SE[, i]
    pred <- Predict[, i]
    
    # Calculate the log-likelihood for this target variable
    log_lik <- sum(dnorm(x = obs, mean = pred, sd = sd, log = TRUE))
    
    # Add a penalty term for extreme deviations
    penalty <- sum((obs - pred)^2)  # You can adjust the penalty term as needed
    
    log_likelihood <- log_likelihood + log_lik - penalty
  }
  
  return(log_likelihood)

}

############################################
## Function that takes a matrix with predictions for each age and sums by 5-years bands

sum_age_predictions <- function(m.data, England.pop) {
  # Create a vector of age groups
  age_groups <- seq(1, nrow(m.data)-1, by = 5)
  
  # Create an empty matrix to store the summarized predictions
  summary_mat <- matrix(0, nrow = length(age_groups), ncol = ncol(m.data))
  
  # Loop through age groups and sum the predictions
  for(n in 1:ncol(m.data)){
  for (i in seq_along(age_groups)) {
    start_age <- age_groups[i]
    end_age <- start_age + 4
    
    summary_mat[i, n] <- mean(m.data[start_age:end_age, n]) #colSums
  }
  }
  # combine old ages to fit the other data (avaiable only for 90+)
  summary_mat_old <- colSums(summary_mat[13:14,])
  summary_mat <- rbind(summary_mat[1:12,], summary_mat_old)
  
  # scale to the England pop size
  summary_mat[,1:7] <- summary_mat[,1:7]*England.pop[1, "Males"]*1000
  summary_mat[,8:14] <- summary_mat[,8:14]*England.pop[1, "Females"]*1000
  rownames(summary_mat)<- NULL
  
  # Return the summary matrix
  return(summary_mat)
}

#######################################

f.GOF.MHA.calc <- function(Targets, SE, Predict){
  
  v.GOF <- rep(0, ncol(Targets)) #create a vector to collect the outputs
  names(v.GOF) <-colnames(Targets)
  
  for(n in 1:ncol(Targets)){
    distribution <- dnorm(x = Targets[,n], mean = Predict[,n], sd = SE[,n], log = T)
    distribution[which(!is.finite(distribution))] <-0 # replace with zero infinite and undefined numbers
    v.GOF[n]<- sum(distribution)  }
  
  # Add more weight to incidence total and mortality
  v.GOF[c(5:6,11:12)]=v.GOF[c(5:6,11:12)]*20
  v.GOF[c(1:2,7:9,14)]=v.GOF[c(1:2,7:9,14)]*10
  
  
  v.GOF
}
###########################################################
# Prior function. A very weak non-influential prior assuming that the uniform distribution in the range of 20% from the fitted value by a random LH approach


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

##############################################################
## Updated prior function in the second calibration attempt

f.prior = function(fitted_params){ #uniform distribution for priors


# Weak priors

    prior.weak <- dunif(fitted_params["P.LGtoHGBC",1], min=fitted_params["P.LGtoHGBC",1]*0.9, max=fitted_params["P.LGtoHGBC",1]*1.1, log = T)
 
# Strong priors
 penalty = -10000000
 
    count.conditions = 
      (fitted_params["P.sympt.diag_LGBC",1] >0.2)*1+
      (fitted_params["P.sympt.diag_LGBC",1] > fitted_params["P.sympt.diag_St1",1])*1+
      (fitted_params["P.sympt.diag_St1",1] > fitted_params["P.sympt.diag_St2",1])*1 +
      (fitted_params["P.sympt.diag_St2",1] > fitted_params["P.sympt.diag_St3",1])*1+
      (fitted_params["P.sympt.diag_St3",1] > fitted_params["P.sympt.diag_St4",1])*1+
      (fitted_params["P.sympt.diag_St1",1] >0.2)*1+
      (fitted_params["P.sympt.diag_St2",1] <0.05)*1+ (fitted_params["P.sympt.diag_St2",1] >0.4)*1+
      (fitted_params["P.sympt.diag_St3",1] <0.1)*1+ (fitted_params["P.sympt.diag_St3",1] >0.7)*1+
      (fitted_params["P.sympt.diag_St4",1] <0.3)*1+ (fitted_params["P.sympt.diag_St3",1] >0.9)*1+
      (fitted_params["P.onset_age",1]<1)*1+ (fitted_params["P.onset_age",1]>1.15)*1+
      (fitted_params["RR.onset_sex",1]<1)*1+ (fitted_params["RR.onset_sex",1]>5)*1+
      (fitted_params["P.sympt.diag_Age",1]< 0.9)*1+
      (fitted_params["shape.t.StI.StII",1] >1)*1+
      (fitted_params["shape.t.StII.StIII",1]>1)*1+
      (fitted_params["shape.t.StIII.StIV",1]>1)*1
    
    prior.strong <- count.conditions*penalty
    
 #prior.strong = if(fitted_params["P.onset",1] <0|
   #fitted_params["P.sympt.diag_LGBC",1] >0.2|
  # fitted_params["P.ungiag.dead",1] < -0.01|
  # fitted_params["P.sympt.diag_LGBC",1] > fitted_params["P.sympt.diag_St1",1]|
  # fitted_params["P.sympt.diag_St1",1] > fitted_params["P.sympt.diag_St2",1]|
  # fitted_params["P.sympt.diag_St2",1] > fitted_params["P.sympt.diag_St3",1]|
  # fitted_params["P.sympt.diag_St3",1] > fitted_params["P.sympt.diag_St4",1]|
  # fitted_params["P.sympt.diag_St1",1] <0 | fitted_params["P.sympt.diag_St1",1] >0.2|
  # fitted_params["P.sympt.diag_St2",1] <0.05 | fitted_params["P.sympt.diag_St2",1] >0.4|
  # fitted_params["P.sympt.diag_St3",1] <0.1 | fitted_params["P.sympt.diag_St3",1] >0.7|
  # fitted_params["P.sympt.diag_St4",1] <0.3 | fitted_params["P.sympt.diag_St3",1] >0.9|
  # fitted_params["P.onset_age",1]<1| fitted_params["P.onset_age",1]>2|
  # fitted_params["RR.onset_sex",1]<1| fitted_params["RR.onset_sex",1]>10|
 #  fitted_params["P.sympt.diag_Age_HGBC",1]< 0.9 |fitted_params["P.sympt.diag_Age_LGBC",1]< 0.9|
  # fitted_params["P.sympt.diag_Age_LGBC",1] > fitted_params["P.sympt.diag_Age_HGBC",1]|
  # fitted_params["P.sympt.diag_Age_HGBC",1] > 1 | fitted_params["P.sympt.diag_Age_LGBC",1] >1|
   #fitted_params["shape.t.StI.StII",1] > fitted_params["shape.t.StII.StIII",1]|
   #fitted_params["shape.t.StII.StIII",1] > fitted_params["shape.t.StIII.StIV",1]|
  # fitted_params["shape.t.StI.StII",1] <0 | fitted_params["shape.t.StI.StII",1] >1|
  # fitted_params["shape.t.StII.StIII",1] <0  |fitted_params["shape.t.StII.StIII",1]>1|
  # fitted_params["shape.t.StIII.StIV",1] <0 | fitted_params["shape.t.StIII.StIV",1]>1){ (-5000000)} else{(0)} #replace with  -Inf
 
return(sum(prior.weak, prior.strong))
}


#############################################################

posterior = function(param){
  return (likelihood(param) + f.prior(param))
}


# The Distributions proposed:
# 4 = Uniform for parameters on age and sex onset and shape for time since onset
# 5 - Normal for C.age.80.undiag.mort
# 7 = Dirichlet for probabilities
# 8 = Truncated normal


####################################################################################
f.proposal.param = function(param, MHA.step){
  
  # Sample for parameters with the uniform distribution

  #param[ ,1] <- rnorm(nrow(param), mean = param[ ,"Mean"], sd= step*param[ ,"Mean"])
  
   param[c(1,2,5:14), 1] <- truncnorm::rtruncnorm(length(param[c(1,2,5:14), 1]), a = 0, b = 1,
                           mean = param[c(1,2,5:14), "Mean"], sd = MHA.step * param[c(1,2,5:14), "Mean"])
  
   param["P.ungiag.dead", "Mean"] <- truncnorm::rtruncnorm(1, a=1, b=-0.01, mean = param["P.ungiag.dead", "Mean"], sd = MHA.step * param["P.ungiag.dead", "Mean"])
    
   
   param[c(3,4) ,"Mean"] <- rnorm(2, mean = param[c(3,4),"Mean"], sd= MHA.step*param[c(3,4),"Mean"])
   

  return(param)
  
}

#############################################
# NOT USED FUNCTIONS
###################################
# Set a strong prior setting all param values between 0 and 1

log_prior_strong <- function(param){
  if((param<0) || (param>1)){  # 
    return(-Inf)}
  else{
    return(0)}
}