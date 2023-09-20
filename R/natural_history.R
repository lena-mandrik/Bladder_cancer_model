### Functions for implementing natural history transitions ###

#' @details
#' This function calculates individualised transition probs for Oc mortality, and BC onset. 
#' @params
#' pop: population matrix containing individual level attributes
#' m.Diag: matrix containing diagnosis status for CRC
#' m.State: matrix containing health states
#' @return matrix of individualised transition probabilities for each transition

f.calc.indiv.TPs <- function(pop, m.Diag, disease){
  
  # Update OCM with all the factors that impact it
  # Other cause mortality transition 
  TP.OC <- OC_mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
  
  # Adjust to reflect an impact of smoking on other cause mortality. Assume that RR for other cause mortality is the same as for all cause mortality 
  TP.OC <- TP.OC*RR.All.Death.no_smoke #RR.All.Death.no_smoke is calibrated 
  TP.OC <- TP.OC*pop[,"no_smoke"] + TP.OC*RR.All.Death.past_smoke*pop[,"past_smoke"] + TP.OC*RR.All.Death.current_smoke*pop[,"current_smoke"]
  
 
  # Retrieve p of cancer onset based on individual factors
  TP.BCLG <- TP.BCHG <- TP.KC <- rep(0, nsample)
  
  if(disease=="bladder" | disease=="bladder_kidney"){
    TP.BCLG <- pop[,"p.BC.i"]*e.BC$P.onset_low.risk
    TP.BCHG <- pop[,"p.BC.i"]*(1-e.BC$P.onset_low.risk)
  }
  
  # Add individual TP for LG to HG cancer based on whether a patient has been diagnosed or not 
  # LG cancers currently don't progress to HG cancers if they are detected unless the period of surveillance passed (2y) in which case they are considered as recurrent LG progressed to HG
  TP.LGtoHGBC <- as.matrix(rep(0, nsample), ncol=1)
  
  if(disease=="bladder" | disease=="bladder_kidney"){
  TP.LGtoHGBC[,1] <- replace(TP.LGtoHGBC[,1], m.Diag[ ,"LG_diag"]==0, e.BC$P.LGtoHGBC) # set probability for those who hasn't been diagnosed yet
  TP.LGtoHGBC[,1] <- replace(TP.LGtoHGBC[,1], m.Diag[ ,"LG_diag"]==1 & (pop[,"age"]-m.Diag[ ,"LG_age_diag"])>1, e.BC$P.LGtoHGBC*e.BC$P.Recurrence.LR) # set probability of recurrence and progression to HG for those who were diagnosed, after 1 year of surveillance when the p is assumed to be zero
  }
  
  TP <- cbind(TP.OC, TP.BCLG, TP.BCHG, TP.LGtoHGBC, pop[,"p.KC.i"])
  colnames(TP) <- c("TP.OC", "TP.BCLG", "TP.BCHG", "TP.LGtoHGBC", "TP.KC")
  
  TP[pop[,"age"] >= 100, "TP.OC"] <- 1
  TP[pop[,"age"] >= 100, c("TP.BCLG", "TP.BCHG", "TP.LGtoHGBC", "TP.KC")] <- 0
  
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
f.calc.Cmort.TP <-function(pop, m.Diag, .env){
  
# Update the transitions for bladder cancer mortality 
TP.C.1.mort <- .env$C1.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.C.1.mort <- TP.C.1.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP.C.2.mort <- .env$C2.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.C.2.mort <- TP.C.2.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP.C.3.mort <- .env$C3.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.C.3.mort <- TP.C.3.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP.C.4.mort <- .env$C4.mort[paste(pop[,"sex"],pop[,"age"], sep = ""), ]
TP.C.4.mort <- TP.C.4.mort[cbind(seq_along(m.Diag[, "yr_diag"]+1), (m.Diag[, "yr_diag"]+1))]

TP_mort <- cbind(TP.C.1.mort, TP.C.2.mort, TP.C.3.mort, TP.C.4.mort)

limit.age <- pop[,"age"] >= 100
TP_mort[limit.age, c("TP.C.1.mort", "TP.C.2.mort", "TP.C.3.mort", "TP.C.4.mort")] <- 0

TP_mort
}

##########################################################################################################################
 #' @details
 #' This function returns whether each individual died from cancer in a sepcific model cycle. 
 #' @params
 #' pop: population matrix containing individual level attributes
 #' m.Diag: matrix containing diagnosis status for CRC
 #' .env: environment for the parameters
 #' m.State: matrix containing health states
 #' @return matrix of individualised values (0 or 1) for each person
 #' 
f.return.C.death <-function(pop, m.Diag, m.State, .env, m.Rand, t){
 
# Calculate individualised transition probs 
TP_C.mort <- f.calc.Cmort.TP(pop, m.Diag, .env)

new_C1_death <- 1*((m.Rand[ ,"Death_C", t] < TP_C.mort[,"TP.C.1.mort"]) & m.State[ ,"St1_HG"]>0)
new_C2_death <- 1*((m.Rand[ ,"Death_C", t] < TP_C.mort[,"TP.C.2.mort"]) & m.State[ ,"St2_HG"]>0)
new_C3_death <- 1*((m.Rand[ ,"Death_C", t] < TP_C.mort[,"TP.C.3.mort"]) & m.State[ ,"St3_HG"]>0)
new_C4_death <- 1*((m.Rand[ ,"Death_C", t] < TP_C.mort[,"TP.C.4.mort"]) & m.State[ ,"St4_HG"]>0)

# TP.ungiag.dead - should be betwen -0.001 and -0.05
# Update the mortality with a possibility for some BC at stage 4 to be undiagnosed death 
died_undiagnosed <- 1*(m.Rand[ ,"Death_C_undiag", t] < (1-exp(.env$P.ungiag.dead * ((pop[, "age"] - 80) * (1*(pop[, "age"] > 80)))))
                       & m.Diag[ ,"yr_diag"]==1
                       & m.Diag[ ,"Sympt_diag"]==1
                       & new_C4_death[]==1)

new_C4_death[died_undiagnosed[] ==1] <-0    

C_death_all <- rowSums(cbind(new_C1_death, new_C2_death, new_C3_death, new_C4_death))

outputs <- list(died_undiagnosed=died_undiagnosed,
                C_death_all =C_death_all)

return(outputs)

}
##########################################################################################################################

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
  
  a.p.it[, 1,]  <- matrix(nrow = nsample, ncol = n.s, c(1 - TP[, "TP.BCLG"]- TP[, "TP.BCHG"]- TP[, "TP.KC"] - TP[, "TP.OC"], TP[, "TP.BCLG"], TP[, "TP.BCHG"], TP[, "TP.KC"], TP[, "TP.OC"]))
  # Transition probabilities in the state No cancer (STATE 1)
  
  a.p.it[, 2,]  <- matrix(nrow = nsample, ncol = n.s, c(rep(0, nsample), 1 - TP[, "TP.LGtoHGBC"] - TP[, "TP.OC"], TP[, "TP.LGtoHGBC"], rep(0, nsample), TP[, "TP.OC"]))
  # Transition probabilities in the state low grade bladder cancer (STATE 2)
  
  a.p.it[, 3,]   <- matrix(nrow = nsample, ncol = n.s, c(rep(0,nsample*2), 1 - TP[, "TP.OC"], rep(0, nsample), TP[, "TP.OC"]))     
  # Transition probabilities in the state high grade bladder cancer (sTATE 3)
  
  a.p.it[, 4,]  <- matrix(nrow = nsample, ncol = n.s, c(rep(0,nsample*3), 1 - TP[, "TP.OC"], TP[, "TP.OC"]))
  # Transition probabilities in the state Kidney cancer (STATE 4)
  
  a.p.it[, 5,]   <- matrix(nrow = nsample, ncol = n.s,  c(rep(0,4),1), byrow = TRUE)
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
  m.Diag[, "yr_onset"][m.Diag[, "yr_onset"] >=1] <- m.Diag[, "yr_onset"][m.Diag[, "yr_onset"] >=1] +1 #Update the year of onset if the cancer developed the previous years
  
  # Add new individuals with HGBC
  new_HG <-  1*(m.M[, t+1] ==3 & m.M[, t] != 3)
  # Add new individuals with LG onset
  new_LG <- m.M[ ,t+1]==2 & m.M[ ,t]==1
  # Add new individuals with KC onset
  new_KC <- m.M[ ,t+1]==4 & m.M[ ,t]==1
  
  # Mark in m.Diag all persons with HGBC, LGBC, KC (independently on diagnosis)
  m.Diag[ ,"HG_state"][which(new_HG>0)] <- 1 # <- replace(m.Diag[ ,"HG_state"], new_HG >0, 1)
  m.Diag[, "LG_state"][which(new_LG>0)] <- 1
  m.Diag[, "KC_state"][which(new_KC>0)] <- 1
  
  # year of onset for each cancer
  m.Diag[, "yr_onset"] <- m.Diag[, "yr_onset"] + new_HG + new_LG + new_KC
  
  # Record age for those individuals who just had BC or KC onset
  m.Diag[, "age_onset"] <- m.Diag[, "age_onset"] + (pop[, "age"] * new_HG)+ (pop[, "age"] * new_KC)   
  m.Diag[, "LG_age_onset"] <- m.Diag[, "LG_age_onset"] + (pop[, "age"] * new_LG) 
  
  return(m.Diag)
}
  

#######################
# Function for allocating HG BC states
############################

f.C.stages <- function(m.M_8s, m.C.T.to.Stage, m.M, m.Diag, t, anal.disease){
  
  # Update m.M_8s matrix with the states for those who have HG cancer
  m.M_8s[, t+1] <- m.M[, t+1]
  
  if(anal.disease =="bladder"){
    m.M_8s[, t+1][which(m.M[, t+1]==4)] <- 1 #Replace with 1 (No bladder cancer) for all with kidney cancer
    c.state= "HG_state"
    n.state=3
  } else {
    m.M_8s[, t+1][which(m.M[, t+1]==2 | m.M[, t+1]==3)] <- 1 #Replace with 1 (NO kidney cancer) for all with bladder cancer
    c.state= "KC_state"
    #for kidney cancer move stage 1 to state 3
    m.M_8s[, t+1][which(m.M[, t+1]==4)] <- 3
    n.state=4
  }
  
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.M_8s[, t] ==8, 8) #Replace with BC death those who died with BC before this cycle
  
  # Replace with the relevant stages for those who were not diagnosed but has cancer
  
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], (m.M[, t+1]==n.state & m.Diag[,c.state]==1 & m.Diag[, "yr_diag"] ==0 & (ceiling(m.C.T.to.Stage[ ,"T.onsetToStage2"])==m.Diag[,"yr_onset"]|m.M_8s[, t]==4)), 4) #Replace for stage 2
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], (m.M[, t+1]==n.state & m.Diag[,c.state]==1 & m.Diag[, "yr_diag"] ==0 & (ceiling(m.C.T.to.Stage[ ,"T.onsetToStage3"])==m.Diag[,"yr_onset"]|m.M_8s[, t]==6)), 6) #Replace for stage 3
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], (m.M[, t+1]==n.state & m.Diag[,c.state]==1 & m.Diag[, "yr_diag"] ==0 & (ceiling(m.C.T.to.Stage[ ,"T.onsetToStage4"])==m.Diag[,"yr_onset"]|m.M_8s[, t]==7)), 7) #Replace for stage 4
  
  # replace with the stage for those who were diagnosed assuming that they don't progress
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "Stage_diag"]==1 & m.Diag[,c.state]==1 & m.M[, t+1] != 5, 3) #replace the stage at diagnosis for those who were diagnosed
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "Stage_diag"]==2 & m.Diag[,c.state]==1 & m.M[, t+1] != 5, 4) #replace the stage at diagnosis for those who were diagnosed
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "Stage_diag"]==3 & m.Diag[,c.state]==1 & m.M[, t+1] != 5, 6) #replace the stage at diagnosis for those who were diagnosed

  return(m.M_8s)
}

