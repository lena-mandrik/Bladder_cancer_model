
#' @details
#' This function sets parameters 
#' @params
#' p.set: Parameter set selected
#' @return all parameters into the global environment
#' 
set_parameters <- function(p.set) {
  
  for(i in 1:length(Param.names)){
    
    vector.Param <- Param_sets[p.set, i]
    names(vector.Param) <- c(Param.names[i])
    assign(names(vector.Param), value =vector.Param, envir = .GlobalEnv)
    
  }
}


#' @details
#' This function generates parameter sets for PSA from parameter distributions 
#' @params
#' Params: Input parameter file containing mean value and distribution information for each model parameter.
#' N_sets: The number of PSA parameter sets to generate
#' @return sets of parameters for PSA. The first row contains mean parameter values for deterministic analysis.


generate_parameters <- function(Params, N_sets){
  
  #Set up matrix for parameter sets
  Param_sets <- matrix(0, ncol = N_sets, nrow = length(Params[,1]))
  rownames(Param_sets) <- rownames(Params)
  Param_sets[,1] <- Params[,1]
  
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
    
    
  }
  
  Param_sets <- t(Param_sets)
  
  ##Test that mean of sample values is roughly similar to mean value - all values should be close to 1
  #Compare_param_means <- Param_sets[1,] / c(apply(Param_sets,2,mean))
  
  Param_sets
  
}

### Functions for processing input data ###

#' @details
#' This function calculates individual bladder cancer risk based on the risk factors 
#' @params
#' pop: population matrix containing individual level attributes
#' Param_sets: File containing sampled parameter sets
#' @return Vector of individual cancer risks
#' 
#Calculate mean population attributes

pop <- population

smoke <- sum(pop[, "current_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])
psmoke <- sum(pop[, "past_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])
occupation <- sum(pop[, "occupation"] * pop[, "weighting"])/sum(pop[, "weighting"])
AsianM <- sum(pop[which(pop[, "ethnic"] ==3 & pop[,"sex"] ==1), "weighting"])/sum(pop[, "weighting"])
AsianF <- sum(pop[which(pop[, "ethnic"] ==3 & pop[,"sex"] ==0), "weighting"])/sum(pop[, "weighting"])
BlackM <- sum(pop[which(pop[, "ethnic"] ==2 & pop[,"sex"] ==1), "weighting"])/sum(pop[, "weighting"])
BlackF <- sum(pop[which(pop[, "ethnic"] ==2 & pop[,"sex"] ==0), "weighting"])/sum(pop[, "weighting"])
OtherM <- sum(pop[which(pop[, "ethnic"] >=4 & pop[,"sex"] ==1), "weighting"])/sum(pop[, "weighting"])
OtherF <- sum(pop[which(pop[, "ethnic"] >=4 & pop[,"sex"] ==0), "weighting"])/sum(pop[, "weighting"])


#Calculate individualised CRC risk by multiplying risks for each attribute
#Note that this equation ensures that weighted product of individual relative risks ==1
#Can't however be sure that weighted product of combined RRs ==1 (not easy to test)

risk <- (RR.current_smoke^(pop[, "current_smoke"] - smoke))*
  (RR.past_smoke ^ (pop[, "past_smoke"] - psmoke))*
  (RR.manufacture ^ (pop[, "occupation"] - occupation))

prod(risk)
no_smoke <- rep(0,nrow(pop))
no_smoke <- replace(no_smoke, pop[, "past_smoke"] !=1 & pop[, "current_smoke"] !=1, 1)
mean(pop[, "current_smoke"]) + mean(pop[, "past_smoke"]) +mean(no_smoke)

pop[, "weighting"]*risk

risk <- (RR_CurrSmoke ^ (pop[, "smo"] - smoke)) *
  (RR_PastSmoke ^ (pop[, "pastsmoke"] - psmoke)) *
  replace((RR_BMI25_M ^ (popA[, "BMI25M"] - BMI25M)), pop[, "sex"] ==0, 1) *
  replace((RR_BMI25_F ^ (popA[, "BMI25F"] - BMI25F)), pop[, "sex"] ==1, 1) *
  replace((RR_BMI30_M ^ (popA[, "BMI30M"] - BMI30M)), pop[, "sex"] ==0, 1) *
  replace((RR_BMI30_F ^ (popA[, "BMI30F"] - BMI30F)), pop[, "sex"] ==1, 1) *
  (RR_LightAlc ^ (popA[, "alcL"] - alcL)) *
  (RR_ModAlc ^ (popA[, "alcM"] - alcM)) *
  (RR_HeavyAlc ^ (popA[, "alcH"] - alcH)) *
  (RR_HighPA ^ (popA[, "PA"] - PA)) *
  (RR_FH ^ (pop[, "FH_real"] - FH)) *
  (RR_Asian_M ^ (popA[, "AsianM"] - AsianM)) *
  (RR_Asian_F ^ (popA[, "AsianF"] - AsianF)) *
  (RR_Black_M ^ (popA[, "BlackM"] - BlackM)) *
  (RR_Black_F ^ (popA[, "BlackF"] - BlackF)) *
  (RR_Other_M ^ (popA[, "OtherM"] - OtherM)) *
  (RR_Other_F ^ (popA[, "OtherF"] - OtherF)) #  *
# (RR_Fibre ^ (pop[, "fibre"] - fibre)) *
# (RR_Meat ^ (pop[, "meat"] - meat)) *

#Make a new temporary matrix where pop attributes match bladder cancer risk categories
popA <- matrix(0, nrow = length(pop[,1]), ncol = 3)

colnames(popA) <- c("BMI25M", "BMI25F", "BMI30M", "BMI30F", "alcL", "alcM", "alcH", "PA", "AsianM", "AsianF", "BlackM", "BlackF", "OtherM", "OtherF")
popA[,"BMI25M"] <- replace(popA[,"BMI25M"], pop[, "bmi"] >=25 & pop[,"bmi"] <30 & pop[,"sex"] ==1, 1)
popA[,"BMI25F"] <- replace(popA[,"BMI25F"], pop[, "bmi"] >=25 & pop[,"bmi"] <30 & pop[,"sex"] ==0, 1)
popA[,"BMI30M"] <- replace(popA[,"BMI30M"], pop[, "bmi"] >=30 & pop[,"sex"] ==1, 1)
popA[,"BMI30F"] <- replace(popA[,"BMI30F"], pop[, "bmi"] >=30 & pop[,"sex"] ==0, 1)
popA[,"alcL"] <- replace(popA[,"alcL"], pop[, "units"] >0 & pop[,"units"] <10.9375, 1)
popA[,"alcM"] <- replace(popA[,"alcM"], pop[, "units"] >10.9375 & pop[,"units"] <43.75, 1)
popA[,"alcH"] <- replace(popA[,"alcH"], pop[, "units"] >43.75, 1)
popA[,"PA"] <- replace(popA[,"PA"], pop[, "PA"] >=600, 1)
popA[,"AsianM"] <- replace(popA[,"AsianM"], pop[, "ethnic"] ==3 & pop[,"sex"] ==1, 1)
popA[,"AsianF"] <- replace(popA[,"AsianF"], pop[, "ethnic"] ==3 & pop[,"sex"] ==0, 1)
popA[,"BlackM"] <- replace(popA[,"BlackM"], pop[, "ethnic"] ==2 & pop[,"sex"] ==1, 1)
popA[,"BlackF"] <- replace(popA[,"BlackF"], pop[, "ethnic"] ==2 & pop[,"sex"] ==0, 1)
popA[,"OtherM"] <- replace(popA[,"OtherM"], pop[, "ethnic"] >=4 & pop[,"sex"] ==1, 1)
popA[,"OtherF"] <- replace(popA[,"OtherF"], pop[, "ethnic"] >=4 & pop[,"sex"] ==0, 1)
