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



 
 #' @details
 #' This function selects the set of transition probabilities that are relevant for each individual given their current health state (probabilities at time t)
 #' @params
 #' m.M: vector containing current health state for each individual
 #' TP: matrix of individualised transition probabilities for each transition
 #' @return matrix of transition probabilities for each individual
 #' 
 #Probs <- function(m.M, TP) {    

 # M_it:    health state occupied by individual i at cycle t (character variable) 
 

 a.p.it <- array(data = NA, dim = c(n.i, n.s, n.s)) # create array of state transition probabilities   
 
 # update m.p.it with the appropriate probabilities (for now there is no probability to die from BC except you are in Clinical state)
 
 a.p.it[, 1,]  <- matrix(nrow = n.i, ncol = n.s, c(1 - TP[, "TP.BC_onset"] - TP[, "TP.OC"] - TP[, "TP.NormToA"], TP[, "TP.NormToLR"], rep(0, n.i), TP[, "TP.NormToA"], rep(0,4*n.i), TP[, "TP.OC"]))
 # Transition probabilities in the state No cancer (STATE 1)
 
 a.p.it[, 2,]  <- matrix(nrow = n.i, ncol = n.s, c(rep(0, n.i), 1 - TP[, "TP.LRtoHR"] - TP[, "TP.OC"], TP[, "TP.LRtoHR"], rep(0, 5*n.i), TP[, "TP.OC"]))
 # Transition probabilities in the state low risk adenoma (STATE 2)
 
 a.p.it[, 3,]   <- matrix(nrow = n.i, ncol = n.s, c(rep(0,n.i*2), 1 - TP[, "TP.HRtoCRCA"] - TP[, "TP.OC"], TP[, "TP.HRtoCRCA"], rep(0, 4*n.i), TP[, "TP.OC"]))  # transition probabilities when sick   
 # Transition probabilities in the state high risk adenoma (sTATE 3)
 
 a.p.it[, 4,]   <- matrix(nrow = n.i, ncol = n.s, c(rep(0,n.i*3), 1  - TP[, "TP.AtoB"] - TP[, "TP.CRC.A.mort"] - TP[, "TP.OC"], TP[, "TP.AtoB"], rep(0,n.i*2), TP[, "TP.CRC.A.mort"], TP[, "TP.OC"]))
 # Transition probabilities in the state proximal preclinical cancer phase I (or A by Duke) (STATE 4)
 
 a.p.it[, 5,]   <- matrix(nrow = n.i, ncol = n.s, c(rep(0,n.i*4), 1  - TP[, "TP.BtoC"] - TP[, "TP.CRC.B.mort"] - TP[, "TP.OC"], TP[, "TP.BtoC"], rep(0,n.i), TP[, "TP.CRC.B.mort"], TP[, "TP.OC"]))
 # Transition probabilities in the state proximal preclinical cancer phase II (or B by Duke) (STATE 5)
 
 a.p.it[, 6,]   <- matrix(nrow = n.i, ncol = n.s, c(rep(0,n.i*5), 1 - TP[, "TP.CtoD"] - TP[, "TP.CRC.C.mort"] - TP[, "TP.OC"], TP[, "TP.CtoD"], TP[, "TP.CRC.C.mort"], TP[, "TP.OC"]))
 # Transition probabilities in the state proximal preclinical cancer phase III (or C by Duke) (STATE 6)
 
 a.p.it[, 7,]   <- matrix(nrow = n.i, ncol = n.s, c(rep(0,n.i*6), 1 - TP[, "TP.CRC.D.mort"] - TP[, "TP.OC"], TP[, "TP.CRC.D.mort"], TP[, "TP.OC"]))
 # Transition probabilities in the state proximal preclinical cancer phase IV (or D by Duke) (STATE 7)
 
 a.p.it[, 8,]   <- matrix(nrow = n.i, ncol = n.s, c(rep(0,7), 1, 0), byrow = TRUE)
 # Transition probabilities in the state death from CRC (STATE 8)
 
 a.p.it[, 9,]   <- matrix(nrow = n.i, ncol = n.s,  c(rep(0,8),1), byrow = TRUE)

# Add to the TP sampling of the time to the next stage and allocation of the stages for everyone who has cancer onset


Mean.t.StI.StII  
Mean.t.StII.StIII 
Mean.t.StIII.StIV

shape.t.StI.StII <- shape.t.StII.StIII <- shape.t.StIII.StIV <- 1.001 

#' @details
#' This function uses the mean and the shape to calculate the scale and sample from the Weibull distribution 
#' 
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
colnames(m.BC.Stage) <-c("T.onsetToStage2", "T.onsetToStage3", "T.onsetToStage4")



m.BC.T.to.Stage <- apply(m.BC.T.to.Stage, 1, FUN =f.stage, Mean.t.StI.StII=Mean.t.StI.StII, 
                         shape.t.StI.StII=shape.t.StI.StII, 
                         Mean.t.StII.StIII =Mean.t.StII.StIII, 
                         shape.t.StII.StIII =shape.t.StII.StIII, 
                         Mean.t.StIII.StIV = Mean.t.StIII.StIV,
                         shape.t.StIII.StIV = shape.t.StIII.StIV)

f.stage(Mean.t.StI.StII, shape.t.StI.StII, Mean.t.StII.StIII, shape.t.StII.StIII, 
        Mean.t.StIII.StIV, shape.t.StIII.StIV)
m.State
m.Diag[ ,"yr_onset"] ==trunc(T.onsetToStage2)