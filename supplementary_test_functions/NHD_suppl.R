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