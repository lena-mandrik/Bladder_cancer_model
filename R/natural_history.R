### Functions for implementing natural history transitions ###

#' @details
#' This function calculates individualised transition probs for Oc mortality, and BC onset. 
#' @params
#' pop: population matrix containing individual level attributes
#' m.Diag: matrix containing diagnosis status for CRC
#' m.State: matrix containing health states
#' @return matrix of individualised transition probabilities for each transition

calc.indiv.TPs <- function(pop, m.Diag){
  
  
  # Other cause mortality transition 
  TP.OC <- OC_mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
  
  # Adjust to reflect an impact of smoking on other cause mortality. Assume that RR for other cause mortality is the same as for all cause mortality 
  TP.OC <- TP.OC*RR.All.Death.no_smoke #RR.All.Death.no_smoke is calibrated 
  TP.OC <- TP.OC*pop[,"no_smoke"] + TP.OC*RR.All.Death.past_smoke*pop[,"past_smoke"] + TP.OC*RR.All.Death.current_smoke*pop[,"current_smoke"]
  
  TP.BC_onset <- pop[,"p.BC.i"]
  
  TP.BCLG <- TP.BC_onset*P.onset_low.risk
  TP.BCHG <- TP.BC_onset*(1-P.onset_low.risk)
  
  # Add individual TP for LG to HG cancer based on whether a patient has been diagnosed or not 
  # LG cancers currently don't progress to HG cancers if they are detected unless the period of surveillance passed (2y) in which case they are considered as recurrent LG progressed to HG
  TP.LGtoHGBC <- as.matrix(rep(0, n.i), ncol=1) 
  TP.LGtoHGBC[,1] <- replace(TP.LGtoHGBC[,1], m.Diag[ ,"LG_diag"]==0, P.LGtoHGBC) # set probability for those who hasn't been diagnosed yet
  TP.LGtoHGBC[,1] <- replace(TP.LGtoHGBC[,1], m.Diag[ ,"LG_diag"]==1 & (pop[,"age"]-m.Diag[ ,"LG_age_diag"])>1, P.LGtoHGBC*P.Recurrence.LR) # set probability of recurrence and progression to HG for those who were diagnosed, after 1 year of surveillance when the p is assumed to be zero
  
  TP <- cbind(TP.OC, TP.BCLG, TP.BCHG, TP.LGtoHGBC)
  colnames(TP) <- c("TP.OC", "TP.BCLG", "TP.BCHG", "TP.LGtoHGBC")
  
  limit.age <- pop[,"age"] >= 100
  TP[limit.age, c("TP.OC", "TP.BCLG", "TP.BCHG", "TP.LGtoHGBC")] <- c(1,0,0,0)

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
  
  a.p.it[, 1,]  <- matrix(nrow = n.i, ncol = n.s, c(1 - TP[, "TP.BCLG"]- TP[, "TP.BCHG"] - TP[, "TP.OC"], TP[, "TP.BCLG"], TP[, "TP.BCHG"], TP[, "TP.OC"]))
  # Transition probabilities in the state No cancer (STATE 1)
  
  a.p.it[, 2,]  <- matrix(nrow = n.i, ncol = n.s, c(rep(0, n.i), 1 - TP[, "TP.LGtoHGBC"] - TP[, "TP.OC"], TP[, "TP.LGtoHGBC"], TP[, "TP.OC"]))
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

#' @details
#' This function updates the matrix with 8 states based on the matrix with 4 states and whether the patient got diagnosed or not
#' It reassigns one of the 4 stages of HG cancer 
#' @params
#' m.M: matrix containing current 4 health state for each individual
#' m.M_8s: matrix containing current 8 health state for each individual
#' @return updated m.M_8s matrix
#' 
f.HG.stage <- function(m.M, m.M_8s, m.BC.T.to.Stage){
  
  m.M_8s[, t+1] <- m.M[, t+1]
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.M_8s[, t] ==8, 8) #Replace with BC death those who died with BC before this cycle
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], ((m.M[, t+1]==3 & m.BC.T.to.Stage[ ,"T.onsetToStage2"]==m.Diag[,"yr_onset"])| (m.M[, t+1]==3 & m.M_8s[, t]==5)) & m.Diag[, "yr_diag"] ==0, 5) #Replace for stage 2
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], ((m.M[, t+1]==3 & m.BC.T.to.Stage[ ,"T.onsetToStage3"]==m.Diag[,"yr_onset"])| (m.M[, t+1]==3 & m.M_8s[, t]==6)) & m.Diag[, "yr_diag"] ==0, 6) #Replace for stage 3
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], ((m.M[, t+1]==3 & m.BC.T.to.Stage[ ,"T.onsetToStage4"]==m.Diag[,"yr_onset"])| (m.M[, t+1]==3 & m.M_8s[, t]==7)) & m.Diag[, "yr_diag"] ==0, 7) #Replace for stage 4
  
  # replace with the stage for those who were diagnosed assuming that they don't progress
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "stage_diag"]==2 & m.M[, t+1] != 4, 5) #replace the stage at diagnosis for those who were diagnosed
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "stage_diag"]==3 & m.M[, t+1] != 4, 6) #replace the stage at diagnosis for those who were diagnosed
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "stage_diag"]==4 & m.M[, t+1] != 4, 7) #replace the stage at diagnosis for those who were diagnosed
  
  m.M_8s
}


#' @details
#' This function updates the time of onset for those persons who had an onset of BC
#' @params
#' m.M: matrix containing current 4 health state for each individual
#' m.Diag: matrix with the characterisics of the diagnosed state
#' @return updated m.M_8s matrix
#' 
f.BC.onset <- function(m.Diag, m.M, pop, t){
  # Record the characteristics of onset for HG cancer
  m.Diag[, "new_onset"] <-0
  m.Diag[, "yr_onset"][m.Diag[, "yr_onset"] >=1] <- m.Diag[, "yr_onset"][m.Diag[, "yr_onset"] >=1] +1 #Update the year of onset if the cancer developed the previous years
  m.Diag[, "new_onset"] <- replace(m.Diag[,"new_onset"], m.M[, t+1] ==3 & m.M[, t] != 3, 1)
  m.Diag[, "yr_onset"] <- m.Diag[, "yr_onset"] + m.Diag[, "new_onset"]
  
  #Mark those individuals who just had BC onset
  m.Diag[, "age_onset"] <- m.Diag[, "age_onset"] + (pop[, "age"] * m.Diag[, "new_onset"])  
  
  # Mark in m.Diag all persons with BC (independently on diagnosis)
  m.Diag[ ,"BC_state"] <- replace(m.Diag[ ,"BC_state"], m.Diag[ ,"yr_onset"] >0, 1)
  
  m.Diag
}


#' @details
#' This function defines who dies from BC
#' @params
#' pop: matrix with population characteristics
#' m.Diag: matrix with the characterisics of the diagnosed state
#' @return updated BC_death_all matrix with 1 if a BC happend at the specific stage
#' 
f.BC.stage.death <- function(pop, m.Diag){
  
  TP_BC.mort <- calc.BCmort.TP(pop, m.Diag)
  
  new_BC1_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.1.mort"]) & m.State[ ,"St1_HG"]>0)
  new_BC2_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.2.mort"]) & m.State[ ,"St2_HG"]>0)
  new_BC3_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.3.mort"]) & m.State[ ,"St3_HG"]>0)
  new_BC4_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.4.mort"]) & m.State[ ,"St4_HG"]>0)
  
  BC_death_all <- rowSums(cbind(new_BC1_death, new_BC2_death, new_BC3_death, new_BC4_death))
  
  BC_death_all
}

#' @details
#' This function records whether a person had a recurrent LG (if they have LGBC) and the number of recurrences 
#' @params
#' P.Recurrence.LR: annual probability of recurrence 
#' m.Diag: matrix with the characterisics of the diagnosed state
#' @return updated m.Diag matrix
#' 
#
f.recurrence.LGBC <- function(m.Diag, n.i, t, P.Recurrence.LR){
  
  TP.BCLG_recurrence <- rep(0, n.i)
  TP.BCLG_recurrence <- replace(TP.BCLG_recurrence, m.Diag[ ,"LG_BC_diag"]==1, P.Recurrence.LR)
  m.Diag[, "LG_BC_n"] <- m.Diag[, "LG_BC_n"] + 1*(m.Rand[ ,"BCLG_recurrence", t] < TP.BCLG_recurrence)
  
  m.Diag
}

 