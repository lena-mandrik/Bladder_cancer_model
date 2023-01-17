### Functions for processing input characteristics for the population and updates them when the functions are called ###

#' @details
#' This function calculates individual probability of onset of BC at time t based on relative risks and probability of onset
#' @params
#' pop: population matrix containing individual level attributes
#' @return Updated pop matrix with individual risk and probabilities of cancer onset in a specific year
#'
#'


f.risk.calc <- function(pop) {
  
  mean.p_smoke <- sum(pop[, "current_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])
  mean.p_psmoke <- sum(pop[, "past_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])
  mean.p_occupation <- sum(pop[, "occupation"] * pop[, "weighting"])/sum(pop[, "weighting"])
  mean.p_sex <- sum(pop[, "sex"] * pop[, "weighting"])/sum(pop[, "weighting"])
  
  #Calculate individual BC risk by multiplying risks for each attribute
  pop[, "risk"] <- (RR.current_smoke^(pop[, "current_smoke"] - mean.p_smoke))*
    (RR.past_smoke ^ (pop[, "past_smoke"] - mean.p_psmoke))*
    (RR.manufacture ^ (pop[, "occupation"] - mean.p_occupation))*(RR.onset_sex^(pop[, "sex"] - mean.p_sex))
  
  # Onset of the bladder cancer at time t by age
  pop[, "p.BC.i"] <-(P.onset*P.onset_age^(pop[, "age"] -30))*pop[, "risk"]
  
  #check
  #no_smoke <- rep(0,nrow(pop))
  #no_smoke <- replace(no_smoke, pop[, "past_smoke"] !=1 & pop[, "current_smoke"] !=1, 1)
  #mean(pop[, "current_smoke"]) + mean(pop[, "past_smoke"]) +mean(no_smoke)
  
  #for (variable in ls()) {
  # assign(variable, get(variable), envir = .GlobalEnv)
  #}
  
  pop
  
}

#' @details
#' This function defines whether a person in HSE remains a smoker based on an annual probability to quit smoking
#' @params
#' pop: population matrix containing individual level attributes
#' @return Updated pop matrix with updated 
#' 

f.smoke.change <- function(pop, rand.quit) {
  # Update smoking status, considering the proportion of smokers who quite smoking during one year
 quit_smoke <- 1*((rand.quit < P.quit.smoke) & pop[,"current_smoke"]==1)
 pop[,"current_smoke"] <- replace(pop[,"current_smoke"],quit_smoke ==1, 0)
 pop[,"past_smoke"] <- replace(pop[,"past_smoke"], pop[,"no_smoke"]==0 & pop[,"past_smoke"]==0 & quit_smoke ==1, 1)
 pop
}