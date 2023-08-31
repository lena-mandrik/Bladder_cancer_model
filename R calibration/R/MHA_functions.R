##  This are the functions used in the Metropolis Hasting algorithm (MHA.R file)

######################################################
# Likelihood function
likelihood = function(Calibr_parameters){
  
  source("R calibration\\Main_model_calibration.R") # run the model
  
   Predict <- f.calibr.output(results_no_screen) # extract the outputs
   output <- cbind(Predict$rate_outcomes_m, Predict$rate_outcomes_f)
  
  #sum by 5 years
  
  #Calculate the predictions by 5-year
  outputs_by_5 <- matrix(0, ncol=(ncol(Targets)), nrow=nrow(Targets))
  
  # get the counter for matrix of data
  n.out=1
  
  for(i.row in 1:(nrow(Targets)-1)){
    for(i.col in 1:ncol(Targets)){
      
      outputs_by_5[i.row,i.col] = mean(output[n.out:(n.out+4),i.col]) #
    }
    n.out=n.out+5
    
  }
  
  #name rows and columns
  colnames(outputs_by_5) <- colnames(Predict)
  rownames(outputs_by_5) <- rownames(Targets)
  

  if(target_stat=="counts"){
   # Predict <- Predict[1:61, ]
    outputs_by_5[ ,1:7] <-outputs_by_5[ ,1:7]*England.pop[,"Males"]*1000
    outputs_by_5[ ,8:14] <-outputs_by_5[ ,8:14]*England.pop[,"Females"]*1000    }
  
  # Exclude mortality from calibration targets
  #if(target_stat=="counts"){
 # Predict <- Predict[1:12, c(1:6,8:13)]
  #m.GOF <- m.GOF[ ,c(1:6,8:13)]
 # SE <- SE[1:12,c(1:6,8:13)]
 # Targets <- Targets[1:12, c(1:6,8:13)]
 # }
  
  likelihood <- f.GOF.MHA.calc(Targets, SE, outputs_by_5)  # log likelihood instead of sum of squared errors
  
  sum_likelihoods = sum(likelihood)
  
  return(sum_likelihoods)

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
  
  # Add more weight to incidence total
  v.GOF[c("incidence.male","incidence.fem")]=v.GOF[c("incidence.male","incidence.fem")]*10
  v.GOF[c("BC.death.male","BC.death.female")]=v.GOF[c("BC.death.male","BC.death.female")]*30
  
  
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

  #  prior.weak <- dunif(fitted_params["P.LGtoHGBC",1], min=fitted_params["P.LGtoHGBC",1]*0.9, max=fitted_params["P.LGtoHGBC",1]*1.1, log = T)
 
# Strong priors
 penalty = -1000
 
    count.conditions = 
      (fitted_params["P.sympt.diag_LGBC",1] >0.2)*1+
      (fitted_params["P.ungiag.dead",1] < -0.01)*1+
      (fitted_params["P.onset_low.risk",1] < 0.6)*1+
      (fitted_params["P.sympt.diag_LGBC",1] > fitted_params["P.sympt.diag_St1",1])*1+
      (fitted_params["P.sympt.diag_St1",1] > fitted_params["P.sympt.diag_St2",1])*1 +
      (fitted_params["P.sympt.diag_St2",1] > fitted_params["P.sympt.diag_St3",1])*1+
      (fitted_params["P.sympt.diag_St3",1] > fitted_params["P.sympt.diag_St4",1])*1+
      (fitted_params["P.sympt.diag_St1",1] >0.2)*1+
      (fitted_params["P.sympt.diag_St2",1] <0.05)*1+ (fitted_params["P.sympt.diag_St2",1] >0.4)*1+
      (fitted_params["P.sympt.diag_St3",1] <0.1)*1+ (fitted_params["P.sympt.diag_St3",1] >0.7)*1+
      (fitted_params["P.sympt.diag_St4",1] <0.3)*1+ (fitted_params["P.sympt.diag_St3",1] >0.9)*1+
      (fitted_params["P.onset_age",1]<1)*1+ (fitted_params["P.onset_age",1]>1.16)*1+
      (fitted_params["RR.onset_sex",1]<1)*1+ (fitted_params["RR.onset_sex",1]>5)*1+
      (fitted_params["P.sympt.diag_Age",1]< 0.9)*1+
      (fitted_params["shape.t.StI.StII",1] >1)*1+
      (fitted_params["shape.t.StII.StIII",1]>1)*1+
      (fitted_params["shape.t.StIII.StIV",1]>1)*1
    
    #Add priors for the mean
    mean.prior1 = (fitted_params["Mean.t.StII.StIII",1] -fitted_params["Mean.t.StI.StII",1])/fitted_params["Mean.t.StI.StII",1]
    mean.prior2 = (fitted_params["Mean.t.StIII.StIV",1] -fitted_params["Mean.t.StII.StIII",1])/fitted_params["Mean.t.StII.StIII",1]
    condition.mean = (mean.prior1 >0.5 | mean.prior2 >0.5)*1
    prior.mean = if(condition.mean==1){-500*mean.prior1+-500*mean.prior2}else{0}
    
    prior <- count.conditions*penalty + prior.mean
    
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
 
return(sum(prior))
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
   
   param["Mean.t.StI.StII", "Mean"] <- truncnorm::rtruncnorm(1, a = 1, b = 7,
                                                  mean = param["Mean.t.StI.StII", "Mean"], sd = MHA.step * param["Mean.t.StI.StII", "Mean"])
   param["Mean.t.StII.StIII", "Mean"] <- truncnorm::rtruncnorm(1, a = 1, b = 5,
                                                             mean = param["Mean.t.StII.StIII", "Mean"], sd = MHA.step * param["Mean.t.StII.StIII", "Mean"])
   param["Mean.t.StIII.StIV", "Mean"] <- truncnorm::rtruncnorm(1, a = 0, b = 4,
                                                             mean = param["Mean.t.StIII.StIV", "Mean"], sd = MHA.step * param["Mean.t.StIII.StIV", "Mean"])
   
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