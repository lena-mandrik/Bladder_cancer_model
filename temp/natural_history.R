### Functions for implementing natural history transitions ###

#' @details
#' This function calculates individualised transition probs for each transition. 
#' @params
#' pop: population matrix containing individual level attributes
#' m.Diag: matrix containing diagnosis status for CRC
#' m.State: matrix containing health states
#' @return matrix of individualised transition probabilities for each transition
#' 
 calc.indiv.TPs <- function(pop, m.Diag, m.State){

# Re-calculate the individual BC risk
pop[, "risk_score"] <- (RR.current_smoke^(pop[, "current_smoke"] - mean.p_smoke))*
  (RR.past_smoke ^ (pop[, "past_smoke"] - mean.p_psmoke))*
  (RR.manufacture ^ (pop[, "occupation"] - mean.p_occupation))*(P.onset_sex^(pop[, "sex"] - mean.p_sex))

# Other cause mortality transition [remember to add the RR for smoking effect]
TP.OC <- OC_mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]

# Adjust to reflect an impact of smoking on other cause mortality. Assume that RR for other cause mortality is the same as for all cause mortality 
TP.OC <- TP.OC*RR.All.Death.no_smoke #RR.All.Death.no_smoke is calibrated 
TP.OC <- TP.OC*pop[,"no_smoke"] + TP.OC*RR.All.Death.past_smoke*pop[,"past_smoke"] + TP.OC*RR.All.Death.current_smoke*pop[,"current_smoke"]

# Onset of the bladder cancer at time t by age
TP.BC_onset <- (P.onset*P.onset_age^(pop[, "age"] -30))*pop[, "risk_score"]

# Probability to become symptomatic patient
# m.Diag[ ,"BC_state"] - should be 1 for each patient with BC in any state.
# Returns an annual probability to be diagnosed by time of cancer onset if a person has cancer

TP.Sympt.diag <- m.Diag[ ,"BC_state"]*P.sympt.diag_A_HGBC*P.sympt.diag_B_HGBC^m.Diag[ ,"yr_onset"]

# Update for those who are older than 80 with the decrement in a symptomatic presentation rate
age_sympt_decrement <- P.sympt.diag_Age80_HGBC^(pop[, "age"]- 80)
age_sympt_decrement <- replace(age_sympt_decrement, age_sympt_decrement>1,1)
TP.Sympt.diag <- TP.Sympt.diag * age_sympt_decrement

# Update the transitions for bladder cancer mortality 

TP.BC.1.mort <- BC.1.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.1.mort <- TP.BC.1.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP.BC.2.mort <- BC.2.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.2.mort <- TP.BC.2.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP.BC.3.mort <- BC.3.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.3.mort <- TP.BC.3.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP.BC.4.mort <- BC.4.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.4.mort <- TP.BC.4.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]


TP <- cbind(TP.BC_onset, TP.Sympt.diag, TP.OC, TP.BC.1.mort, TP.BC.2.mort, TP.BC.3.mort, TP.BC.4.mort)

limit.age <- pop[,"age"] >= 100
TP[limit.age, c("TP.BC_onset", "TP.Sympt.diag", "TP.OC", "TP.BC.1.mort", "TP.BC.2.mort", "TP.BC.3.mort", "TP.BC.4.mort")] <- 0
TP[limit.age, "TP.OC"] <- 1

TP

}

Mean.t.StI.StII 
Mean.t.StII.StIII 
Mean.t.StIII.StIV

shape.t.StI.StII <- shape.t.StII.StIII <- shape.t.StIII.StIV <- 1.001 



# Add to the TP sampling of the time to the next stage and allocation of the stages for everyone who has cancer onset

w_mean =4
shape =1.5

f.stage <- function(w_mean, shape){
  scale = w_mean/gamma(1+1/shape)
  
  time_to_state <- rweibull(1, shape, scale)
  
  time_to_state
  
}

f.stage(w_mean, shape.t.StI.StII)
