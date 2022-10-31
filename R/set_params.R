########################################################
## Functions to allocate the time to stage at diagnosis for each person in HSE
#' @details
#' This function uses the mean and the shape to calculate the scale and sample from the Weibull distribution 

f.Weibull.sample <- function(w_mean, shape){
  scale = w_mean/gamma(1+1/shape)
  
  time_to_state <- rweibull(1, shape, scale)
  
  time_to_state
  
}

#' @details
#' This function returns a time from onset to the Stage 2,3, and 4 for each person in the model
#' @params
#' Mean.t. - mean time from each stage till the next one. Based on a fixed input from a qualitative study
#' shape.t. - a calibrated shape for each Weibull distribution for each parameter
#' @return matrix of transition probabilities for each individual

f.stage <- function(Mean.t.StI.StII, shape.t.StI.StII, Mean.t.StII.StIII, shape.t.StII.StIII, 
                    Mean.t.StIII.StIV, shape.t.StIII.StIV){
  
  T.onsetToStage2 <- f.Weibull.sample(Mean.t.StI.StII, shape.t.StI.StII)
  T.Stage3 <- f.Weibull.sample(Mean.t.StII.StIII, shape.t.StII.StIII)
  T.Stage4 <- f.Weibull.sample(Mean.t.StIII.StIV, shape.t.StIII.StIV)
  
  T.onsetToStage3 <- T.onsetToStage2+T.Stage3
  T.onsetToStage4 <- T.onsetToStage3+T.Stage4
  v.Stages <- c(T.onsetToStage2, T.onsetToStage3, T.onsetToStage4)
  
  v.Stages
  
}

#create a matrix for all people in the dataset with the time they get each stage if they get invasive cancer
m.BC.T.to.Stage <- matrix(nrow = n.i, ncol = 3)
colnames(m.BC.T.to.Stage) <-c("T.onsetToStage2", "T.onsetToStage3", "T.onsetToStage4")
rownames(m.BC.T.to.Stage) <- 1:n.i

#' @details
#' This code assigns the time till next stage for each person in HSE probabilistically
#' @params 
#' m.BC.T.to.Stage - matrix containing the time to the stage for each person in HSE

f.stage.assign <- function(m.BC.T.to.Stage){
  
  for(i in 1:n.i){
    m.BC.T.to.Stage[i,] <- f.stage(Mean.t.StI.StII, shape.t.StI.StII, Mean.t.StII.StIII, shape.t.StII.StIII, Mean.t.StIII.StIV, shape.t.StIII.StIV)
  }
  m.BC.T.to.Stage <- round(m.BC.T.to.Stage)
  m.BC.T.to.Stage
}

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
#' @return Updated pop matrix with individual risk and probabilities of cancer onset in a specific year
#'
#'


f.risk.calc <- function(pop) {
  
#Calculate mean population attributes.
# Sex associated risk is the only calibrated parameter
  
mean.p_smoke <- sum(pop[, "current_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])
mean.p_psmoke <- sum(pop[, "past_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])
mean.p_occupation <- sum(pop[, "occupation"] * pop[, "weighting"])/sum(pop[, "weighting"])
mean.p_sex <- sum(pop[, "sex"] * pop[, "weighting"])/sum(pop[, "weighting"])

#Calculate individual BC risk by multiplying risks for each attribute
pop[, "risk_score"] <- (RR.current_smoke^(pop[, "current_smoke"] - mean.p_smoke))*
  (RR.past_smoke ^ (pop[, "past_smoke"] - mean.p_psmoke))*
  (RR.manufacture ^ (pop[, "occupation"] - mean.p_occupation))*(P.onset_sex^(pop[, "sex"] - mean.p_sex))

#TP.BC_onset <- (P.onset*P.onset_age^(pop[, "age"] -30))*pop[, "risk_score"]

#check
#no_smoke <- rep(0,nrow(pop))
#no_smoke <- replace(no_smoke, pop[, "past_smoke"] !=1 & pop[, "current_smoke"] !=1, 1)
#mean(pop[, "current_smoke"]) + mean(pop[, "past_smoke"]) +mean(no_smoke)

for (variable in ls()) {
  assign(variable, get(variable), envir = .GlobalEnv)
}

pop

}

#' @details
#' This function sets up an array of random numbers
#' @params
#' @return an array of random numbers for each event, person and time cycle
generate_random <- function() {
  events <- c("PROBS", "SYMPT", "Death_BC", "Screen_UPTK","Diag_UPTK")
  array(runif(n.i * length(events) * n.t), dim = c(n.i, length(events), n.t), dimnames = list(NULL, events, NULL))
} 


#' @details
#' This function calculates individual other cause mortality based on the risk factors 
#' @params
#' pop: population
#' 
 #Set up random number array for each individual
m.Rand <- generate_random()


