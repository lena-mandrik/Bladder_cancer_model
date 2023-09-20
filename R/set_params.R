
#################################################################################################################################
# Create the disease-specific environments to save the future parameters
#################################################################################################################################
if(disease == "kidney" | disease == "bladder_kidney"){e.KC <- new.env()}
if(disease == "bladder" | disease == "bladder_kidney"){e.BC <- new.env()}


#################################################################################################################################
# Weibull sampling forprogression of undiagnosed cancer
#################################################################################################################################
## Functions to allocate the time to the next stage for each person in HSE
#' @details
#' This function uses the mean and the shape to calculate the scale and sample from the Weibull distribution 


f.Weibull.sample <- function(w_mean, shape.p, nsample){
  
  scale.p = w_mean/gamma(1+1/shape.p)
  
  time_to_state <- rweibull(nsample, shape.p, scale.p)
  
  return(time_to_state)
  
}


#################################################################################################################################
# Individual sampling of the speed of cancer progression based on Weibull distribution
#################################################################################################################################

#' @details
#' This function returns a time from onset to the Stage 2,3, and 4 for each person in the model
#' @params
#' Mean.t. - mean time from each stage till the next one. Based on a fixed input from a qualitative study
#' shape.t. - a calibrated shape for each Weibull distribution for each parameter
#' @return matrix of transition probabilities for each individual

f.stage <- function(.env, nsample){
  
  
  Mean.t.StI.StII =.env$Mean.t.StI.StII
  shape.t.StI.StII=.env$shape.t.StI.StII
  Mean.t.StII.StIII=.env$Mean.t.StII.StIII
  shape.t.StII.StIII=.env$shape.t.StII.StIII
  Mean.t.StIII.StIV=.env$Mean.t.StIII.StIV
  shape.t.StIII.StIV=.env$shape.t.StIII.StIV

  T.onsetToStage2 <- f.Weibull.sample(Mean.t.StI.StII, shape.t.StI.StII, nsample)
  T.Stage3 <- f.Weibull.sample(Mean.t.StII.StIII, shape.t.StII.StIII, nsample)
  T.Stage4 <- f.Weibull.sample(Mean.t.StIII.StIV, shape.t.StIII.StIV, nsample)
  
  # Adjust T.Stage3 and T.Stage4 based on the value of T.onsetToStage2 (i.e. that quicker growing cancers grow quicker at all stages)
  T.Stage3 <- T.Stage3 * (T.onsetToStage2 / mean(T.onsetToStage2))
  T.Stage4 <- T.Stage4 * (T.onsetToStage2 / mean(T.onsetToStage2))
  
  T.onsetToStage3 <- T.onsetToStage2+T.Stage3
  T.onsetToStage4 <- T.onsetToStage3+T.Stage4
  v.Stages <- cbind(T.onsetToStage2, T.onsetToStage3, T.onsetToStage4) #
  
  v.Stages
  
}

#################################################################################################################################
# Cancer survival
#################################################################################################################################
#' @details
#' This function sets the survival for either bladder or kidney cancers 
#' @params
#' OC_BC_mort: a txt table with survival data
#' @return all parameters into the disease specific environment
#' 
#' 

f.C.mort <- function(list_of_file_names, set.envir=.GlobalEnv, disease_name){
  
  
  list_of_files <- lapply(list_of_file_names, read.table, sep = "\t", header =T) #Add the files to the list
  
  names(list_of_files) <- tools::file_path_sans_ext(basename(list_of_file_names))  # Save the names
  
  list2env(list_of_files,envir=set.envir) # Save to the relevant environment 
  
  # Proceeding cancer mortality
  mort1 <- matrix(0, ncol = 1, nrow = 142)
  mort2 <- matrix(0, ncol = 60, nrow = 142)
  C1.mort <- cbind(mort1, as.matrix(set.envir$S1), mort2)
  C2.mort <- cbind(mort1, as.matrix(set.envir$S2), mort2)
  C3.mort <- cbind(mort1, as.matrix(set.envir$S3), mort2)
  C4.mort <- cbind(mort1, as.matrix(set.envir$S4), mort2)
  
  colnames(C1.mort) <- colnames(C2.mort) <- colnames(C3.mort) <- colnames(C4.mort) <- c(0:70)
  rownames(C1.mort) <- rownames(C2.mort) <- rownames(C3.mort) <- rownames(C4.mort) <- c(paste("0",c(30:100), sep = ""), paste("1",c(30:100), sep = ""))
  
  # Select the objects for the outputs
  output <- list(C1.mort=C1.mort, C2.mort=C2.mort, C3.mort=C3.mort, C4.mort=C4.mort)
  
  for (variable in names(output)) {
    #new_variable_name <- paste(variable, disease_name, sep = "_")
    assign(variable, output[[variable]], envir = set.envir)
      }
}


#################################################################################################################################
# Setting up disease-specific parameters
#################################################################################################################################
#' @details
#' This function sets parameters 
#' @params
#' p.set: Parameter set selected
#' @return all parameters into the selected environment
#' 
#' 
f.set_parameters <- function(p.set, Param_sets, Param.names, set.envir=.GlobalEnv, disease_name) { #choose "_KC", "_BC"
  
  for(i in 1:26){ #length(Param.names) for NHD parameters
    
    vector.Param <- Param_sets[p.set, i]
    #names(vector.Param) <- paste0(c(Param.names[i]), disease_name)
    names(vector.Param) <- c(Param.names[i])                      
    assign(names(vector.Param), value =vector.Param, envir =set.envir)
  }
  
  # set parameters for the symptomatic diagnosis
  Symp_params <- c(0,
                   Param_sets[p.set, "P.sympt.diag_LGBC"],
                   Param_sets[p.set, "P.sympt.diag_St1"],
                   Param_sets[p.set, "P.sympt.diag_St2"],
                   0,
                   Param_sets[p.set, "P.sympt.diag_St3"],
                   Param_sets[p.set, "P.sympt.diag_St4"],
                   0)
  
  
  #Set the parameters for tests sensitivity and False positive (from specificity) in a format of a value by state to use in matrix multiplication
  test_accuracy <- as.matrix(c((1-Param_sets[p.set, "Spec.dipstick"]), #for no cancer
                     Param_sets[p.set, "Sens.dipstick.LG"], #LG state
                     Param_sets[p.set, "Sens.dipstick.St1"], #stage 1
                     Param_sets[p.set, "Sens.dipstick.St2.4"], #stage 2
                     0, # Death OC
                     rep(Param_sets[p.set, "Sens.dipstick.St2.4"],2), #stages 3-4
                     0), ncol=1) # Death cancer
                     
  #Set the parameters for the diagnostic accuracy of US +cytology (assumed to be the same for all HG cancers)               
  diag1_accuracy <- as.matrix(c((1-Param_sets[p.set, "Spec.joined.diag"]), # for no cancer
                                Param_sets[p.set, "Sens.joined.diag.LG"], # LG state
                                Param_sets[p.set, "Sens.joined.diag.HG"], # stage 1
                                Param_sets[p.set, "Sens.joined.diag.HG"], # stage 2
                                0, # Death OC
                                rep(Param_sets[p.set, "Sens.joined.diag.HG"],2), #stages 3-4
                                0), ncol=1) # Death cancer
  
  #Set the parameters for the diagnostic accuracy of flexible cystoscopy (assumed to be the same for all HG cancers)               
  diag2_accuracy <- as.matrix(c((1-Param_sets[p.set, "Spec.cystoscopy"]), # for no cancer
                               Param_sets[p.set, "Sens.cystoscopy.LG"], # LG state
                               Param_sets[p.set, "Sens.cystoscopy.HG"], # stage 1
                               Param_sets[p.set, "Sens.cystoscopy.HG"], # stage 2
                                0, # Death OC
                                rep(Param_sets[p.set, "Sens.cystoscopy.HG"],2), #stages 3-4
                                0), ncol=1) # Death cancer
  
  rownames(test_accuracy)  <- rownames(diag1_accuracy)  <- rownames(diag2_accuracy)  <- states_long
  colnames(test_accuracy) <- colnames(diag1_accuracy) <-colnames(diag2_accuracy) <-"Sens"
  
  # Set harms: mortality due to TURBT                  
  Mort.TURBT <- Param_sets[p.set, "Mort.TURBT"]
  
  # Set utility age and stage decrements
  Disutility.HG.St1.3 <- Param_sets[p.set, "Disutility.HG.St1.3"]
  Disutility.HG.St4 <- Param_sets[p.set, "Disutility.HG.St4"]
  Disutility.LG <- Param_sets[p.set, "Disutility.LG"]
  
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
  


  for (variable in ls()) {
    assign(variable, get(variable), envir = set.envir)
  }
 }


#################################################################################################################################
# Setting up general parameters identical for both diseases
#################################################################################################################################
#' @details
#' This function sets parameters that are identical for bladder and kidney models
#' @params
#' p.set: Parameter set selected
#' @return all parameters into the Global  environment
#' 
#' 
f.set_gen_parameters <- function(p.set, Param_sets, Param.names) { 
  
  for(i in 1:4){ #length(Param.names) for NHD parameters
    
    vector.Param <- Param_sets[p.set, i]
    #names(vector.Param) <- paste0(c(Param.names[i]), disease_name)
    names(vector.Param) <- c(Param.names[i])                      
    assign(names(vector.Param), value =vector.Param, envir =.GlobalEnv)
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

  Utility.age <- Param_sets[p.set, "Utility.age"]

  # Create Cost matrix summarising the cost parameters for each screening process
  m.Cost.screen <- as.matrix(c(Param_sets[p.set, "Cost.ad.dipstick"],
                             Param_sets[p.set, "Cost.dipstick.positive"],
                             Param_sets[p.set, "Cost.dipstick"]), ncol=1)

  # Set matrix for screening and diagnostic costs
  m.scr.diag.costs <- as.matrix(c(Param_sets[p.set, "Cost.diag.sympt"],
                                  Param_sets[p.set, "Cost.diag.screen1"],
                                  Param_sets[p.set, "Cost.diag.screen2"],
                                  Param_sets[p.set, "Cost.dipstick.invite"],
                                  Param_sets[p.set, "Cost.ad.dipstick"],
                                  Param_sets[p.set, "Cost.dipstick.positive"],
                                  Param_sets[p.set, "Cost.dipstick"]))
  
  DS_names <-rownames(m.Cost.screen) <- c("Invite_DS","Respond_DS", "Positive_DS")
  
  # Create matrix summarising the screening process
  screen_names <- c("Invite_DS","Respond_DS", "Positive_DS",
                    "Respond_diag", "Positive_diag", "Respond_Cyst", 
                    "Diagnostic_Cyst", "Surgery",
                    "Next_Surv", 
                    "Die_Surgery", "FP_BC", "FP_KC", "FN_BC", "FN_KC",
                    "HG", "LG", "KC")
  
  for (variable in ls()) {
    assign(variable, get(variable), envir = .GlobalEnv)
  }

}

#################################################################################################################################
# Generating parameter sets
#################################################################################################################################
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

#################################################################################################################################
# Generating random numbers for each individual
#################################################################################################################################

#' @details
#' This function sets up an array of random numbers
#' @params
#' @return an array of random numbers for each event, person and time cycle
f.generate_random <- function() {
  events <- c("Screen_time", "PROBS", "Smoke_quit", "BCLG_recurrence", "SYMPT", "Death_C", "Respond_DS", "Positive_DS", 
              "Respond_diag", "Positive_diag", "Respond_Cyst", "Positive_Cyst", "Die_Surgery", "Death_C_undiag")
  array(runif(nsample * length(events) * n.t), dim = c(nsample, length(events), n.t), dimnames = list(NULL, events, NULL))
} 
