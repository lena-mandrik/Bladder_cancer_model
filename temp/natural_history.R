### Functions for implementing natural history transitions ###

#' @details
#' This function calculates individualised transition probs for each transition. 
#' @params
#' pop: population matrix containing individual level attributes
#' m.Diag: matrix containing diagnosis status for CRC
#' m.State: matrix containing health states
#' @return matrix of individualised transition probabilities for each transition
#' 
# calc.indiv.TPs <- function(pop, m.Diag, m.State){

pop[, "risk_score"] <- (RR.current_smoke^(pop[, "current_smoke"] - mean.p_smoke))*
  (RR.past_smoke ^ (pop[, "past_smoke"] - mean.p_psmoke))*
  (RR.manufacture ^ (pop[, "occupation"] - mean.p_occupation))*(P.onset_sex^(pop[, "sex"] - mean.p_sex))

# Other cause mortality transition [remember to add the RR for smoking effect]
TP.OC <- OC_mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]

# Onset of the bladder cancer at time t by age
TP.BC_onset <- (P.onset*P.onset_age^(pop[, "age"] -30))*pop[, "risk_score"]


# Probability to become symptomatic patient
# m.Diag[ ,"BC_state"] - should be 1 for each patient with BC in any state.
# Returns an annual probability to be diagnosed by time of cancer onset if a person has cancer

TP.Sympt.diag <- m.Diag[ ,"BC_state"]*P.sympt.diag_A_HGBC*P.sympt.diag_B_HGBC^m.Diag[ ,"yr_onset"]

# update for those who are older than 80 with the decrement in a symptomatic presentation rate
age_sympt_decrement <- P.sympt.diag_Age80_HGBC*(pop[, "age"]- 80)
age_sympt_decrement <- replace(age_sympt_decrement, age_sympt_decrement>0,0)
TP.Sympt.diag <- TP.Sympt.diag + age_sympt_decrement


TP.CRC.A.mort <- CRC.A.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.CRC.A.mort <- TP.CRC.A.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]
TP.CRC.B.mort <- CRC.B.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.CRC.B.mort <- TP.CRC.B.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]
TP.CRC.C.mort <- CRC.C.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.CRC.C.mort <- TP.CRC.C.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]
TP.CRC.D.mort <- CRC.D.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]