########################################################
## Functions to allocate the time to stage at diagnosis for each person in HSE
#' @details
#' This function uses the mean and the shape to calculate the scale and sample from the Weibull distribution 


f.Weibull.sample <- function(w_mean, shape, nsample){
  
  scale = w_mean/gamma(1+1/shape)
  
  time_to_state <- rweibull(nsample, shape, scale)
  
  time_to_state
  
}

#' @details
#' This function returns a time from onset to the Stage 2,3, and 4 for each person in the model
#' @params
#' Mean.t. - mean time from each stage till the next one. Based on a fixed input from a qualitative study
#' shape.t. - a calibrated shape for each Weibull distribution for each parameter
#' @return matrix of transition probabilities for each individual

f.stage <- function(Mean.t.StI.StII, shape.t.StI.StII, Mean.t.StII.StIII, shape.t.StII.StIII, 
                    Mean.t.StIII.StIV, shape.t.StIII.StIV, nsample){
  
  T.onsetToStage2 <- f.Weibull.sample(Mean.t.StI.StII, shape.t.StI.StII, nsample)
  T.Stage3 <- f.Weibull.sample(Mean.t.StII.StIII, shape.t.StII.StIII, nsample)
  T.Stage4 <- f.Weibull.sample(Mean.t.StIII.StIV, shape.t.StIII.StIV, nsample)
  
  T.onsetToStage3 <- T.onsetToStage2+T.Stage3
  T.onsetToStage4 <- T.onsetToStage3+T.Stage4
  #v.Stages <- round(cbind(T.onsetToStage2, T.onsetToStage3, T.onsetToStage4)) #
  v.Stages <- cbind(T.onsetToStage2, T.onsetToStage3, T.onsetToStage4) #
  
  v.Stages
  
}


#' @details
#' This function sets parameters 
#' @params
#' p.set: Parameter set selected
#' @return all parameters into the global environment
#' 
#' 
f.set_parameters <- function(p.set) {
  
  for(i in 1:25){ #length(Param.names) for NHD parameters
    
    vector.Param <- Param_sets[p.set, i]
    names(vector.Param) <- c(Param.names[i])
    
    assign(names(vector.Param), value =vector.Param, envir = .GlobalEnv)
  }
  
  #Set the parameters for Dipstick uptake
  coef_DT_Uptk <- c(Param_sets[p.set, "DT.UPTK.CONS"],
                     Param_sets[p.set, "DT.UPTK.50"],
                     Param_sets[p.set, "DT.UPTK.55"],
                     Param_sets[p.set, "DT.UPTK.65"],
                     Param_sets[p.set, "DT.UPTK.70"],
                     Param_sets[p.set, "DT.UPTK.F"],
                     Param_sets[p.set, "DT.UPTK.NRESP"],
                     Param_sets[p.set, "DT.UPTK.INC"],
                     Param_sets[p.set, "DT.UPTK.IMD2"],
                     Param_sets[p.set, "DT.UPTK.IMD3"],
                     Param_sets[p.set, "DT.UPTK.IMD4"],
                     Param_sets[p.set, "DT.UPTK.IMD5"],
                     Param_sets[p.set, "DT.UPTK.ASIAN"])
  
  #Set the parameters for diagnostic uptake with screen positive
  Diag.UPTK <- Param_sets[p.set, "Diag.UPTK"]
  
  #Set the parameters for tests sensitivity and False positive (from specificity) in a format of a value by state to use in matrix multiplication
  test_accuracy <- as.matrix(c((1-Param_sets[p.set, "Spec.dipstick"]), #for no cancer
                     Param_sets[p.set, "Sens.dipstick.LG"], #LG state
                     Param_sets[p.set, "Sens.dipstick.St1"], #stage 1
                     0, # Death OC
                     rep(Param_sets[p.set, "Sens.dipstick.St2.4"],3), #stages 2-4
                     0), ncol=1) # Death BC
                     
  #Set the parameters for the diagnostic accuracy of US +cytology (assumed to be the same for all HG cancers)               
  diag1_accuracy <- as.matrix(c((1-Param_sets[p.set, "Spec.joined.diag"]), #for no cancer
                                Param_sets[p.set, "Sens.joined.diag.LG"], #LG state
                                Param_sets[p.set, "Sens.joined.diag.HG"], #stage 1
                                0, # Death OC
                                rep(Param_sets[p.set, "Sens.joined.diag.HG"],3), #stages 2-4
                                0), ncol=1) # Death BC
  
  #Set the parameters for the diagnostic accuracy of flexible cystoscopy (assumed to be the same for all HG cancers)               
  diag2_accuracy <- as.matrix(c((1-Param_sets[p.set, "Spec.cystoscopy"]), #for no cancer
                               Param_sets[p.set, "Sens.cystoscopy.LG"], #LG state
                               Param_sets[p.set, "Sens.cystoscopy.HG"], #stage 1
                              0, # Death OC
                             rep(Param_sets[p.set, "Sens.cystoscopy.HG"],3), #stages 2-4
                               0), ncol=1) # Death BC  
  
  rownames(test_accuracy)  <- rownames(diag1_accuracy)  <- rownames(diag2_accuracy)  <- states_long
  colnames(test_accuracy) <- colnames(diag1_accuracy) <-colnames(diag2_accuracy) <-"Sens"
  
  # Set harms: mortality due to TURBT                  
  Mort.TURBT <- Param_sets[p.set, "Mort.TURBT"]
  
  # Set utility age and stage decrements
  Utility.age <- Param_sets[p.set, "Utility.age"]
  Disutility.HG.St1.3 <- Param_sets[p.set, "Disutility.HG.St1.3"]
  Disutility.HG.St4 <- Param_sets[p.set, "Disutility.HG.St4"]
  Disutility.LG <- Param_sets[p.set, "Disutility.LG"]
  
  # Set costs of diagnosis and screening
  Cost.diag.sympt <- Param_sets[p.set, "Cost.diag.sympt"]
  Cost.diag.screen1 <- Param_sets[p.set, "Cost.diag.screen1"]
  Cost.diag.screen2 <- Param_sets[p.set, "Cost.diag.screen2"]
  
  # Create Cost matrix summarising the cost parameters for each screening process
  m.Cost.screen <- as.matrix(c(Param_sets[p.set, "Cost.ad.dipstick"],
                   Param_sets[p.set, "Cost.dipstick.positive"],
                   Param_sets[p.set, "Cost.dipstick"]), ncol=1)
  
  DS_names <-rownames(m.Cost.screen) <- c("Invite_DS","Respond_DS", "Positive_DS")
  
  # Create matrix summarising the screening process
  screen_names <- c("Invite_DS","Respond_DS", "Positive_DS","Respond_diag", "Positive_diag", "Respond_Cyst", "Diagnostic_Cyst", "TURBT",
                    "Next_Surv", 
                    "Die_TURBT", "FP", "FN",
                    "HG", "LG")

  # Set treatment and surveillance costs
  # Create Cost matrices 
  m.Cost.treat <- matrix(0, nrow =n.t+1, ncol =n.s_long)
  rownames(m.Cost.treat) <- c(0:n.t)
  colnames(m.Cost.treat) <- states_long
  
  #For LG
  m.Cost.treat["1", "BC_LG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.LG"]
  m.Cost.treat["2", "BC_LG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.LG"]+Param_sets[p.set, "Cost.treat.Y2"]
  m.Cost.treat["3", "BC_LG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.LG"]+Param_sets[p.set, "Cost.treat.Y3"]
  
  #For Stage 1
  m.Cost.treat["1", "St1_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage1"]
  m.Cost.treat["2", "St1_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage1"]+Param_sets[p.set, "Cost.treat.Y2"]
  m.Cost.treat["3", "St1_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage1"]+Param_sets[p.set, "Cost.treat.Y3"]

  #For Stage 2
  m.Cost.treat["1", "St2_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage2"]
  m.Cost.treat["2", "St2_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage2"]+Param_sets[p.set, "Cost.treat.Y2"]
  m.Cost.treat["3", "St2_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage2"]+Param_sets[p.set, "Cost.treat.Y3"]
  
  #For Stage 3
  m.Cost.treat["1", "St3_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage3"]
  m.Cost.treat["2", "St3_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage3"]+Param_sets[p.set, "Cost.treat.Y2"]
  m.Cost.treat["3", "St3_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage3"]+Param_sets[p.set, "Cost.treat.Y3"]
  
  #For Stage 4
  m.Cost.treat["1", "St4_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage4"]
  m.Cost.treat["2", "St4_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage4"]+Param_sets[p.set, "Cost.treat.Y2"]
  m.Cost.treat["3", "St4_HG"] <- Param_sets[p.set, "Cost.treat.intercept"] + Param_sets[p.set, "Cost.treat.stage4"]+Param_sets[p.set, "Cost.treat.Y3"]
  
  #Add surveillance costs for Y4,5 for HG cancers only
  m.Cost.treat[c("4","5"), c("St1_HG","St2_HG","St3_HG","St4_HG")] <- Param_sets[p.set, "Cost.surv.Y4.5"]
  
  # Set matrix for screening and diagnostic costs
  m.scr.diag.costs <- as.matrix(c(Param_sets[p.set, "Cost.diag.sympt"],
                                Param_sets[p.set, "Cost.diag.screen1"],
                                Param_sets[p.set, "Cost.diag.screen2"],
                                Param_sets[p.set, "Cost.dipstick.invite"],
                                Param_sets[p.set, "Cost.ad.dipstick"],
                                Param_sets[p.set, "Cost.dipstick.positive"],
                                Param_sets[p.set, "Cost.dipstick"]))
                                

  for (variable in ls()) {
    assign(variable, get(variable), envir = .GlobalEnv)
  }
 }


#' @details
#' This function generates parameter sets for PSA from parameter distributions 
#' @params
#' Params: Input parameter file containing mean value and distribution information for each model parameter.
#' N_sets: The number of PSA parameter sets to generate
#' @return sets of parameters for PSA. The first row contains mean parameter values for deterministic analysis.


f.generate_parameters <- function(Params, N_sets){
  
  #Set up matrix for parameter sets
  Param_sets <- matrix(0, ncol = N_sets, nrow = length(Params[,1]))
  rownames(Param_sets) <- rownames(Params)
  Param_sets[,1] <- Params[,"Param3"]
  
  if(N_sets > 1){
    
    #Sample from distributions 1=Beta; 2=Gamma, 3=Log-normal, 4=Uniform, 5=Normal, 6=Constant
    
    Param_sets[which(Params[, "Distribution"] ==1),2:N_sets] <- rbeta(length(Param_sets[which(Params[, "Distribution"] ==1),1])*(N_sets-1), 
                                                                      shape1=Params[which(Params[, "Distribution"] ==1), "Param1"], 
                                                                      shape2=Params[which(Params[, "Distribution"] ==1), "Param2"])
    
    Param_sets[which(Params[, "Distribution"] ==2),2:N_sets] <- rgamma(length(Param_sets[which(Params[, "Distribution"] ==2),1])*(N_sets-1), 
                                                                       shape=Params[which(Params[, "Distribution"] ==2), "Param1"], 
                                                                       scale=Params[which(Params[, "Distribution"] ==2), "Param2"])
    
    Param_sets[which(Params[, "Distribution"] ==3),2:N_sets] <- rlnorm(length(Param_sets[which(Params[, "Distribution"] ==3),1])*(N_sets-1), 
                                                                       meanlog=Params[which(Params[, "Distribution"]==3), "Param1"], 
                                                                       sdlog=Params[which(Params[, "Distribution"] ==3), "Param2"])
    
    Param_sets[which(Params[, "Distribution"] ==4),2:N_sets] <- runif(length(Param_sets[which(Params[, "Distribution"] ==4),1])*(N_sets-1), 
                                                                      min=Params[which(Params[, "Distribution"] ==4), "Param1"], 
                                                                      max=Params[which(Params[, "Distribution"] ==4), "Param2"])
    
    Param_sets[which(Params[, "Distribution"] ==5),2:N_sets] <- rnorm(length(Param_sets[which(Params[, "Distribution"] ==5),1])*(N_sets-1), 
                                                                      mean=Params[which(Params[, "Distribution"] ==5), "Param1"], 
                                                                      sd=Params[which(Params[, "Distribution"] ==5), "Param2"])
    
    Param_sets[which(Params[, "Distribution"] ==6),2:N_sets] <- Param_sets[which(Params[, "Distribution"] ==6),1]
    
    Param_sets[which(Params[, "Distribution"] ==7),2:N_sets] <- truncnorm::rtruncnorm(length(Param_sets[which(Params[, "Distribution"] ==7),1])*(N_sets-1), 
                                                                      a=0, b=1, mean=Params[which(Params[, "Distribution"] ==7), "Param1"], 
                                                                      sd=Params[which(Params[, "Distribution"] ==7), "Param2"])   
    # truncated distribution is used for some calibrated parameters that represents porbabilties. They were restricted at [0,1] interval at calibration  by assigning -Inf priors and are truncated here

  }
  
  Param_sets <- t(Param_sets)
  
  ##Test that mean of sample values is roughly similar to mean value - all values should be close to 1
  #Compare_param_means <- Param_sets[1,] / c(apply(Param_sets,2,mean))
  
  Param_sets
  
}



#' @details
#' This function sets up an array of random numbers
#' @params
#' @return an array of random numbers for each event, person and time cycle
f.generate_random <- function() {
  events <- c("Screen_time", "PROBS", "Smoke_quit", "BCLG_recurrence", "SYMPT_HG", "SYMPT_LG","Death_BC", "Respond_DS", "Positive_DS", 
              "Respond_diag", "Positive_diag", "Respond_Cyst", "Positive_Cyst", "Die_TURBT")
  array(runif(nsample * length(events) * n.t), dim = c(nsample, length(events), n.t), dimnames = list(NULL, events, NULL))
} 
