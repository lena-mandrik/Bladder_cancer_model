### Functions for implementing natural history transitions ###

#' @details
#' This function calculates individualised transition probs for Oc mortality, and BC onset. 
#' @params
#' pop: population matrix containing individual level attributes
#' m.Diag: matrix containing diagnosis status for CRC
#' m.State: matrix containing health states
#' @return matrix of individualised transition probabilities for each transition

f.calc.indiv.TPs <- function(pop, m.Diag, BC.4.mort){
  
  # Update OCM with all the factors that impact it
  
  # Other cause mortality transition 
  TP.OC <- OC_mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
  
  # Adjust to reflect an impact of smoking on other cause mortality. Assume that RR for other cause mortality is the same as for all cause mortality 
  TP.OC <- TP.OC*RR.All.Death.no_smoke #RR.All.Death.no_smoke is calibrated 
  TP.OC <- TP.OC*pop[,"no_smoke"] + TP.OC*RR.All.Death.past_smoke*pop[,"past_smoke"] + TP.OC*RR.All.Death.current_smoke*pop[,"current_smoke"]
  
 
  # Retrieve p of cancer onset based on individual factors
  TP.BC_onset <- pop[,"p.BC.i"]
  TP.BCLG <- TP.BC_onset*P.onset_low.risk
  TP.BCHG <- TP.BC_onset*(1-P.onset_low.risk)
  
  # Add individual TP for LG to HG cancer based on whether a patient has been diagnosed or not 
  # LG cancers currently don't progress to HG cancers if they are detected unless the period of surveillance passed (2y) in which case they are considered as recurrent LG progressed to HG
  TP.LGtoHGBC <- as.matrix(rep(0, nsample), ncol=1) 
  TP.LGtoHGBC[,1] <- replace(TP.LGtoHGBC[,1], m.Diag[ ,"LG_diag"]==0, P.LGtoHGBC) # set probability for those who hasn't been diagnosed yet
  TP.LGtoHGBC[,1] <- replace(TP.LGtoHGBC[,1], m.Diag[ ,"LG_diag"]==1 & (pop[,"age"]-m.Diag[ ,"LG_age_diag"])>1, P.LGtoHGBC*P.Recurrence.LR) # set probability of recurrence and progression to HG for those who were diagnosed, after 1 year of surveillance when the p is assumed to be zero
  
  TP <- cbind(TP.OC, TP.BCLG, TP.BCHG, TP.LGtoHGBC)
  colnames(TP) <- c("TP.OC", "TP.BCLG", "TP.BCHG", "TP.LGtoHGBC")
  
  TP[pop[,"age"] >= 100, "TP.OC"] <- 1
  TP[pop[,"age"] >= 100, c("TP.BCLG", "TP.BCHG", "TP.LGtoHGBC")] <- 0
  
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
f.calc.BCmort.TP <-function(pop, m.Diag){
  
# Update the transitions for bladder cancer mortality 
TP.BC.1.mort <- BC.1.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.1.mort <- TP.BC.1.mort[cbind(seq_along(m.Diag[, "HG_yr_diag"]+1), (m.Diag[, "HG_yr_diag"]+1))]

TP.BC.2.mort <- BC.2.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.2.mort <- TP.BC.2.mort[cbind(seq_along(m.Diag[, "HG_yr_diag"]+1), (m.Diag[, "HG_yr_diag"]+1))]

TP.BC.3.mort <- BC.3.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.3.mort <- TP.BC.3.mort[cbind(seq_along(m.Diag[, "HG_yr_diag"]+1), (m.Diag[, "HG_yr_diag"]+1))]

TP.BC.4.mort <- BC.4.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.BC.4.mort <- TP.BC.4.mort[cbind(seq_along(m.Diag[, "HG_yr_diag"]+1), (m.Diag[, "HG_yr_diag"]+1))]

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
f.Probs <- function(m.M, TP) {    
  
  # M_it:    health state occupied by individual i at cycle t (character variable) 
  
  
  a.p.it <- array(data = NA, dim = c(nsample, n.s, n.s)) # create array of state transition probabilities   
  
  # update m.p.it with the appropriate probabilities (for now there is no probability to die from BC except you are in Clinical state)
  
  a.p.it[, 1,]  <- matrix(nrow = nsample, ncol = n.s, c(1 - TP[, "TP.BCLG"]- TP[, "TP.BCHG"] - TP[, "TP.OC"], TP[, "TP.BCLG"], TP[, "TP.BCHG"], TP[, "TP.OC"]))
  # Transition probabilities in the state No cancer (STATE 1)
  
  a.p.it[, 2,]  <- matrix(nrow = nsample, ncol = n.s, c(rep(0, nsample), 1 - TP[, "TP.LGtoHGBC"] - TP[, "TP.OC"], TP[, "TP.LGtoHGBC"], TP[, "TP.OC"]))
  # Transition probabilities in the state low grade bladder cancer (STATE 2)
  
  a.p.it[, 3,]   <- matrix(nrow = nsample, ncol = n.s, c(rep(0,nsample*2), 1 - TP[, "TP.OC"], TP[, "TP.OC"]))     
  # Transition probabilities in the state high grade bladder cancer (sTATE 3)
  
  a.p.it[, 4,]   <- matrix(nrow = nsample, ncol = n.s,  c(rep(0,3),1), byrow = TRUE)
  # Transition probabilities in the state mortality
  
  # Add to the TP sampling of the time to the next stage and allocation of the stages for everyone who has cancer onset
  
  m.p.it <- matrix(NA, nrow = nsample, ncol = n.s)
  for (i in 1:nsample){
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
f.samplev <- function (m.p, m.Rand, t) {
  lev <- 1:n.s
  U <- apply(m.p, 1, cumsum)
  
  if (any((U[n.s, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  rand <- rep(m.Rand[, "PROBS", t], rep(n.s, nsample))
  new <- lev[1 + colSums(rand > U)]
  
  new
}



#' @details
#' This function records whether a person had a recurrent LG (if they have LGBC) and the number of recurrences 
#' @params
#' P.Recurrence.LR: annual probability of recurrence 
#' m.Diag: matrix with the characterisics of the diagnosed state
#' @return updated m.Diag matrix
#' 
#
f.recurrence.LGBC <- function(m.Diag, nsample, t, P.Recurrence.LR){
  
  TP.BCLG_recurrence <- rep(0, nsample)
  TP.BCLG_recurrence <- replace(TP.BCLG_recurrence, m.Diag[ ,"LG_BC_diag"]==1, P.Recurrence.LR)
  m.Diag[, "LG_BC_n"] <- m.Diag[, "LG_BC_n"] + 1*(m.Rand[ ,"BCLG_recurrence", t] < TP.BCLG_recurrence)
  
  m.Diag
}






###########

f.onset.Diag <- function(m.Diag, m.M, pop, t){
  
  # Record the characteristics of onset for HG cancer
  m.Diag[, "HG_yr_onset"][m.Diag[, "HG_yr_onset"] >=1] <- m.Diag[, "HG_yr_onset"][m.Diag[, "HG_yr_onset"] >=1] +1 #Update the year of onset if the cancer developed the previous years
  new_HG <-  1*(m.M[, t+1] ==3 & m.M[, t] != 3)
  m.Diag[, "HG_yr_onset"] <- m.Diag[, "HG_yr_onset"] + new_HG
  
  # Mark those individuals who just had BC onset
  m.Diag[, "HG_age_onset"] <- m.Diag[, "HG_age_onset"] + (pop[, "age"] * new_HG)  
  
  # Mark in m.Diag all persons with BC (independently on diagnosis)
  m.Diag[ ,"HG_state"] <- replace(m.Diag[ ,"HG_state"], m.Diag[ ,"HG_yr_onset"] >0, 1)
  
  # Mark age of those individuals who just had LG onset
  new_LG <- m.M[ ,t+1]==2 & m.M[ ,t]==1
  
  # Record the characteristics of onset for LG cancer
  m.Diag[, "LG_state"][which(new_LG)] <- 1
  m.Diag[, "LG_age_onset"] <- m.Diag[, "LG_age_onset"] + (pop[, "age"] * new_LG) 
  
  return(m.Diag)
}
  

#######################
f.HG.stages <- function(m.M_8s, m.BC.T.to.Stage, m.M, m.Diag, t){
  
  # Update m.M_8s matrix with the states for those who have HG cancer
  m.M_8s[, t+1] <- m.M[, t+1]
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.M_8s[, t] ==8, 8) #Replace with BC death those who died with BC before this cycle
  
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], ((m.M[, t+1]==3 & ceiling(m.BC.T.to.Stage[ ,"T.onsetToStage2"])==m.Diag[,"HG_yr_onset"])| (m.M[, t+1]==3 & m.M_8s[, t]==5)) & m.Diag[, "HG_yr_diag"] ==0, 5) #Replace for stage 2
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], ((m.M[, t+1]==3 & ceiling(m.BC.T.to.Stage[ ,"T.onsetToStage3"])==m.Diag[,"HG_yr_onset"])| (m.M[, t+1]==3 & m.M_8s[, t]==6)) & m.Diag[, "HG_yr_diag"] ==0, 6) #Replace for stage 3
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], ((m.M[, t+1]==3 & ceiling(m.BC.T.to.Stage[ ,"T.onsetToStage4"])==m.Diag[,"HG_yr_onset"])| (m.M[, t+1]==3 & m.M_8s[, t]==7)) & m.Diag[, "HG_yr_diag"] ==0, 7) #Replace for stage 4
  
  # replace with the stage for those who were diagnosed assuming that they don't progress
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "HG_stage_diag"]==2 & m.M[, t+1] != 4, 5) #replace the stage at diagnosis for those who were diagnosed
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "HG_stage_diag"]==3 & m.M[, t+1] != 4, 6) #replace the stage at diagnosis for those who were diagnosed
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "HG_stage_diag"]==4 & m.M[, t+1] != 4, 7) #replace the stage at diagnosis for those who were diagnosed
  
  return(m.M_8s)
}

##################################
# Function updating the matrix m.BC.T.to.Stage with shorter or longer progression based whether the age is below or above 60

#f.HG.stages <- function(m.BC.T.to.Stage, pop){
 # speed_by_age =1.0001
  
  # Set up benchmark as 60
 # age_diff_from_benchmark = pop[,"age"] - 60
  
 # age_adjustment_factor = ifelse(age_diff_from_benchmark >= 0, speed_by_age*age_diff_from_benchmark, abs(age_diff_from_benchmark) * (1-(speed_by_age-1)))
  
#  m.BC.T.to.Stage= m.BC.T.to.Stage*(pop[ ,"age"]-60)*speed_by_age
#  m.BC.T.to.Stage
#}