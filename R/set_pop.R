### Functions for processing input characteristics for the population and updates them when the functions are called ###

#' @details
#' This function calculates individual probability of onset of BC at time t based on relative risks and probability of onset
#' @params
#' pop: population matrix containing individual level attributes
#' @return Updated pop matrix with individual risk and probabilities of cancer onset in a specific year
#'
#' The possible env are _KC for kidney cancer and _BC for bladder cancer


f.risk.calc <- function(pop, set.envir, disease) {
  
  mean.p_smoke <- sum(pop[, "current_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])
  mean.p_psmoke <- sum(pop[, "past_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])
  mean.p_occupation <- sum(pop[, "occupation"] * pop[, "weighting"])/sum(pop[, "weighting"])
  mean.p_sex <- sum(pop[, "sex"] * pop[, "weighting"])/sum(pop[, "weighting"])
  
  #extract from the environment
  RR.current_smoke=set.envir$RR.current_smoke
  RR.past_smoke=set.envir$RR.past_smoke
  RR.manufacture=set.envir$RR.manufacture
  RR.onset_sex=set.envir$RR.onset_sex
  P.onset=set.envir$P.onset
  P.onset_age=set.envir$P.onset_age
  
  #Calculate individual BC or KC  risk by multiplying risks for each attribute
  
  if(disease=="kidney"){
    pop[, "p.KC.i"] <- ((RR.current_smoke^(pop[, "current_smoke"] - mean.p_smoke))*
                        (RR.past_smoke ^ (pop[, "past_smoke"] - mean.p_psmoke))*
                        (RR.manufacture ^ (pop[, "occupation"] - mean.p_occupation))*(RR.onset_sex^(pop[, "sex"] - mean.p_sex)))*(P.onset*P.onset_age^(pop[, "age"] -30))
    
  } else {
    pop[, "p.BC.i"] <- ((RR.current_smoke^(pop[, "current_smoke"] - mean.p_smoke))*
                          (RR.past_smoke ^ (pop[, "past_smoke"] - mean.p_psmoke))*
                          (RR.manufacture ^ (pop[, "occupation"] - mean.p_occupation))*(RR.onset_sex^(pop[, "sex"] - mean.p_sex)))*(P.onset*P.onset_age^(pop[, "age"] -30))
    
  }
  
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
    return(pop)

}

####################################################

#' @details
#' This function calculates the smoking status of the former smokers based on their age when the age is re-set to earlier
#' @params
#' pop: population matrix containing individual level attributes
#' @return Updated pop matrix with updated 
#' 

f.smoke.initial <- function(pop, m.Rand) {
  
  #define for how many years the smoking status should be re-calculated
  years_dif = pop[,"age_0"]-pop[,"age"]
  
  years_dif_cum=matrix(0, nrow=nrow(pop), ncol=1)
  
  for(i in 1:nrow(pop)){
    years_dif_cum[i,1]  = sum(m.Rand[i,"Smoke_quit", 1:years_dif[i]])
  }
  
  # Update smoking status, considering the proportion of smokers who quite smoking during one year
  status_smoke <- 1*((years_dif_cum[,1] > P.quit.smoke) & pop[,"past_smoke"]==1)
  pop[,"past_smoke"] <- replace(pop[,"past_smoke"],status_smoke ==1, 0)
  pop[,"current_smoke"] <- replace(pop[,"current_smoke"], pop[,"no_smoke"]==0 & pop[,"past_smoke"]==0 & status_smoke ==1, 1)
  
  return(pop)
  
}


########################################################
set_population <- function(population, disease) {
  
  if(cohort ==0){ #if cohort ==0, individuals start in model at true (HSE) age
    
    
# Re-assign the weights if the population of a specific age range is selected 
    weights_England <- as.matrix(population[ ,"weighting"], ncol=1)
    weights_England[,1] <- replace(weights_England[,1], population[ ,"age_0"]> max_age | population[ ,"age_0"]< min_age, 0)

# Sample according to the desired sample size
    samp_idx <- sample(seq_len(nrow(population)), nsample, replace = TRUE, prob=weights_England[,1])
    new_population <- population[samp_idx, ]
    state <- rep(0, nsample)
    p.States <- states_pop[new_population[,"PID"], ]

    for(i in 1:nsample) {
      x <-colnames(p.States)
      state[i] <- as.numeric(sample(x, size=1, replace = T, prob = p.States[i,]))
    }
    
    state <- replace(state, state==9, 5) # assign those with the diagnosed cancer as "dead" (ie exclude from costs /benefits calculation)
    population <- cbind(new_population, state) 
    
    # Exclude population with diagnosed cancer by assigning them to the dead state (8)
    
    if(disease == "bladder_kidney" | disease == "bladder"){population <-f.risk.calc(population, e.BC, "bladder")}
    if(disease == "bladder_kidney" | disease == "kidney"){population <-f.risk.calc(population, e.KC, "kidney")}
     } else{
       if(disease == "bladder_kidney" | disease == "bladder"){population <-f.risk.calc(population, e.BC, "bladder")}
       if(disease == "bladder_kidney" | disease == "kidney"){population <-f.risk.calc(population, e.KC, "kidney")}
    }
  population
    }
