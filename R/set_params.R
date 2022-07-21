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