### Functions for implementing natural history transitions ###

#' @details
#' This function calculates individualised transition probs for Oc mortality, and BC onset. 
#' @params
#' pop: population matrix containing individual level attributes
#' m.Diag: matrix containing diagnosis status for CRC
#' m.State: matrix containing health states
#' @return matrix of individualised transition probabilities for each transition

 calc.indiv.TPs <- function(pop, m.Diag, m.State){

# Re-calculate the individual BC risk
pop[, "risk_score"] <- (RR.current_smoke^(pop[, "current_smoke"] - mean.p_smoke))*
  (RR.past_smoke ^ (pop[, "past_smoke"] - mean.p_psmoke))*
  (RR.manufacture ^ (pop[, "occupation"] - mean.p_occupation))*(P.onset_sex^(pop[, "sex"] - mean.p_sex))

# Other cause mortality transition 
TP.OC <- OC_mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]

# Adjust to reflect an impact of smoking on other cause mortality. Assume that RR for other cause mortality is the same as for all cause mortality 
TP.OC <- TP.OC*RR.All.Death.no_smoke #RR.All.Death.no_smoke is calibrated 
TP.OC <- TP.OC*pop[,"no_smoke"] + TP.OC*RR.All.Death.past_smoke*pop[,"past_smoke"] + TP.OC*RR.All.Death.current_smoke*pop[,"current_smoke"]

# Onset of the bladder cancer at time t by age
TP.BC_onset <- (P.onset*P.onset_age^(pop[, "age"] -30))*pop[, "risk_score"]

TP <- cbind(TP.BC_onset, TP.OC)

limit.age <- pop[,"age"] >= 100
TP[limit.age, c("TP.BC_onset", "TP.OC")] <- 0
TP[limit.age, "TP.OC"] <- 1

TP
 }
 
 #' @details
 #' This function calculates individualised transition probs for BC mortality. 
 #' @params
 #' pop: population matrix containing individual level attributes
 #' m.Diag: matrix containing diagnosis status for CRC
 #' m.State: matrix containing health states
 #' @return matrix of individualised transition probabilities for each transition
 #' 
calc.BCmort.TP <-function(pop, m.Diag){
  
# Update the transitions for bladder cancer mortality 
TP.BC.1.mort <- BC.1.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.1.mort <- TP.BC.1.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP.BC.2.mort <- BC.2.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.2.mort <- TP.BC.2.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP.BC.3.mort <- BC.3.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.3.mort <- TP.BC.3.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP.BC.4.mort <- BC.4.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.4.mort <- TP.BC.4.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP_mort <- cbind(TP.BC.1.mort, TP.BC.2.mort, TP.BC.3.mort, TP.BC.4.mort)

limit.age <- pop[,"age"] >= 100
TP_mort[limit.age, c("TP.BC.1.mort", "TP.BC.2.mort", "TP.BC.3.mort", "TP.BC.4.mort")] <- 0

TP_mort
}

#' @details
#' This function defines who dies during the cycle from BC
#' @params
#' 

#TP_BC_death <- calc.BCmort.TP(pop, m.Diag)

#Determine who is newly died
#new_BC_death_1 <- (m.Rand[ ,"Death_BC", t] < TP_BC_death[,"TP.BC.1.mort"]) & (m.Diag[, "BC_diag"] ==1 & m.Diag[, "stage_diag"] ==1)
#new_BC_death_2 <- (m.Rand[ ,"Death_BC", t] < TP_BC_death[,"TP.BC.2.mort"]) & (m.Diag[, "BC_diag"] ==1 & m.Diag[, "stage_diag"] ==2)
#new_BC_death_3 <- (m.Rand[ ,"Death_BC", t] < TP_BC_death[,"TP.BC.3.mort"]) & (m.Diag[, "BC_diag"] ==1 & m.Diag[, "stage_diag"] ==3)
#new_BC_death_4 <- (m.Rand[ ,"Death_BC", t] < TP_BC_death[,"TP.BC.4.mort"]) & (m.Diag[, "BC_diag"] ==1 & m.Diag[, "stage_diag"] ==4)



 #' @details
 #' This function selects the set of transition probabilities that are relevant for each individual given their current health state (probabilities at time t)
 #' @params
 #' m.M: vector containing current health state for each individual
 #' TP: matrix of individualised transition probabilities for each transition
 #' @return matrix of transition probabilities for each individual
 #' 
Probs <- function(m.M, TP) {    

 # M_it:    health state occupied by individual i at cycle t (character variable) 
 

 a.p.it <- array(data = NA, dim = c(n.i, n.s, n.s)) # create array of state transition probabilities   
 
 # update m.p.it with the appropriate probabilities (for now there is no probability to die from BC except you are in Clinical state)
 
 a.p.it[, 1,]  <- matrix(nrow = n.i, ncol = n.s, c(1 - TP[, "TP.BC_onset"] - TP[, "TP.OC"], TP[, "TP.BC_onset"]*P.onset_low.risk, TP[, "TP.BC_onset"]*(1-P.onset_low.risk), TP[, "TP.OC"]))
 # Transition probabilities in the state No cancer (STATE 1)
 
 a.p.it[, 2,]  <- matrix(nrow = n.i, ncol = n.s, c(rep(0, n.i), 1 - rep(P.LGtoHGBC, n.i) - TP[, "TP.OC"], rep(P.LGtoHGBC, n.i), TP[, "TP.OC"]))
 # Transition probabilities in the state low grade bladder cancer (STATE 2)
 
 a.p.it[, 3,]   <- matrix(nrow = n.i, ncol = n.s, c(rep(0,n.i*2), 1 - TP[, "TP.OC"], TP[, "TP.OC"]))     
 # Transition probabilities in the state high grade bladder cancer (sTATE 3)
 
 a.p.it[, 4,]   <- matrix(nrow = n.i, ncol = n.s,  c(rep(0,3),1), byrow = TRUE)
 # Transition probabilities in the state mortality

# Add to the TP sampling of the time to the next stage and allocation of the stages for everyone who has cancer onset
 
 m.p.it <- matrix(NA, nrow = n.i, ncol = n.s)
 for (i in 1:n.i){
   m.p.it[i,] <- a.p.it[i, m.M[i], ]
 }
 
 m.p.it
 }



#' @details
#' This function implements the state transition and updates health states
#' @params
#' m.p: matrix of transition probabilities for each individual
#' m.Rand: matrix of individual random numbers
#' t: current cycle
#' @return a vector of new health states
samplev <- function (m.p, m.Rand, t) {
  lev <- 1:n.s
  U <- apply(m.p, 1, cumsum)
  
  if (any((U[n.s, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  rand <- rep(m.Rand[, "PROBS", t], rep(n.s, n.i))
  new <- lev[1 + colSums(rand > U)]
  
  new
}



