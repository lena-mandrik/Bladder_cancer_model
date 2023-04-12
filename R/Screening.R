### Functions for screening ###

#' @details
#' This function calculates individualised screening parameters for each person
#' @params
#' pop: a population matrix of individuals, each with a current age
#' m:Screen: an array of screening history
#' m.Diag: a matrix giving diagnostic status for all individuals
#' m.State: a matrix giving current health state for all individuals
#' @return a matrix of individualised screening parameters
#' 

calc.screen.Params <- function(pop, m.Screen, m.Diag, m.State) {
  
    #Calculates dipstick uptake by personal characteristics
    DS_uptake <- 1/(1+exp(-((cbind((m.State[, "DeathBC"] ==0 & m.State[, "DeathOC"] ==0), pop[,"age"] <55, 
                                    pop[,"age"] >=55 & pop[,"age"] <60, pop[,"age"] >=65 & pop[,"age"] <70, pop[,"age"] >=70, 
                                    pop[, "sex"] ==0, (rowSums(m.Screen[, ,"Invite_DS"]) >=1 & rowSums(m.Screen[, ,"Respond_DS"]) ==0),
                                    (rowSums(m.Screen[, ,"Respond_DS"]) >=1), pop[, "imd"] ==2, pop[, "imd"] ==3, pop[, "imd"] ==4, 
                                    pop[, "imd"] ==5, pop[, "ethnic"] ==3) * 1) %*% coef_DT_Uptk)))
    
    #Calculates DT sensitivity or false positives depending upon underlying health state
    DS_diag <- (m.State %*% test_accuracy[, "Sens"]) 
    

  scr.Params <- cbind(DS_uptake, DS_diag)
  
  colnames(scr.Params) <- c("DS_uptake", "DS_diag")
  
  scr.Params
}



#' @details
#' This function runs FIT or gFOBT screening
#' @params
#' m:Screen: an array of screening and surveillance history
#' m.Diag: a matrix giving diagnostic status for all individuals
#' m.State: a matrix of current health states for all individuals
#' m.Rand: a random number array
#' pop: a population matrix of individuals, each with a current age
#' t: current time point
#' scr.Params: updated screening parameters
#' FITage: age at which biennial FIT screening starts
#' @return an updated screening history array
#' 
DS_screen <- function(m.Screen, m.Diag, m.State, m.Rand, pop, t, scr.Params, FITage) {
  
  #determine who is eligible for FIT screening
  FIT_elig <- m.Diag[, "CRC_diag"] ==0 & m.Screen[ ,t , "Next_Surv"] ==0 & 
    m.State[, "DeathCRC"] ==0 & m.State[, "DeathOC"] ==0 &
    (pop[, "age"] >= FITage | pop[, "age"] >= pop[, "risk_age"]) & 
    pop[, "age"] <=74 & m.Screen[ ,t , "Invite_FIT"] ==0
  
  #invite eligible
  m.Screen[ ,t+1 , "Invite_FIT"] <- replace(m.Screen[ ,t+1 , "Invite_FIT"],  FIT_elig ==TRUE, 1)
  
  #decide who takes up invite
  FIT_invite <- (m.Rand[ ,"FIT_UPTK", t] < scr.Params[, "FIT_uptake"]) & FIT_elig
  m.Screen[ ,t+1 , "Respond_FIT"] <- replace(m.Screen[ ,t+1 , "Respond_FIT"],  FIT_invite ==TRUE, 1)
  
  #Include a repeat FIT in inadequates DON'T NEED THIS AS ALREADY IN COSTINGS
  #m.Screen[ ,t+1 , "Respond_FIT"] <- m.Screen[ ,t+1 , "Respond_FIT"] + ((m.Rand[ ,"FIT_INAD", t] < inad_FIT) & FIT_invite) * 1 
  
  #decide who gets a positive result (includes true positives and false positives)
  FIT_Pos <- ((m.Rand[ ,"FIT_POS", t] < scr.Params[, "FIT_prob_pos"]) & FIT_invite)
  m.Screen[ ,t+1 , "Positive_FIT"] <- replace(m.Screen[ ,t+1 , "Positive_FIT"],  FIT_Pos ==TRUE, 1)
  
  #decide who is invited to colonoscopy and who is invited to CTC
  CTC_invite <- (m.Rand[ ,"COL_CTC", t] < scr.Params[, "CTC_%"]) & m.Screen[ ,t+1 , "Positive_FIT"] ==1
  m.Screen[ ,t+1 , "Invite_CTC"] <- replace(m.Screen[ ,t+1 , "Invite_CTC"], CTC_invite ==TRUE, 1)
  Col_invite <- (m.Rand[ ,"COL_CTC", t] >= scr.Params[, "CTC_%"]) & m.Screen[ ,t+1 , "Positive_FIT"] ==1
  m.Screen[ ,t+1 , "Invite_Scrn_Col"] <- replace(m.Screen[ ,t+1 , "Invite_Scrn_Col"], Col_invite ==TRUE, 1)
  
  #Decide who takes up CTC and colonoscopy
  CTC_attend <- (m.Rand[ ,"CTC_UPTK", t] < CTC_Uptk) & CTC_invite
  Col_attend <- (m.Rand[ ,"COL_UPTK", t] < scr.Params[, "Col_FIT_uptake"]) & Col_invite
  
  #diagnose those found through CTC/Colonoscopy screening
  LR <- ((m.Rand[ ,"COL_POS", t] < COL_LR_sens) & Col_attend & m.State[, "LR"] ==1) |
    ((m.Rand[ ,"CTC_POS", t] < CTC_LR_sens) & CTC_attend & m.State[, "LR"] ==1)
  HR <- ((m.Rand[ ,"COL_POS", t] < COL_HR_sens) & Col_attend & m.State[, "HR"] ==1) |
    ((m.Rand[ ,"CTC_POS", t] < CTC_HR_sens) & CTC_attend & m.State[, "HR"] ==1) 
  CRC <- ((m.Rand[ ,"COL_POS", t] < COL_CRC_sens) & Col_attend &
            (m.State[, "CRCA"] ==1 | m.State[, "CRCB"] ==1 | m.State[, "CRCC"] ==1 | m.State[, "CRCD"] ==1)) |
    ((m.Rand[ ,"CTC_POS", t] < CTC_CRC_sens) & CTC_attend &
       (m.State[, "CRCA"] ==1 | m.State[, "CRCB"] ==1 | m.State[, "CRCC"] ==1 | m.State[, "CRCD"] ==1))
  
  m.Screen[ ,t+1 , "LR_Adenoma"] <- replace(m.Screen[ ,t+1 , "LR_Adenoma"], LR ==TRUE, 1)
  m.Screen[ ,t+1 , "HR_Adenoma"] <- replace(m.Screen[ ,t+1 , "HR_Adenoma"], HR ==TRUE, 1)
  m.Screen[ ,t+1 , "CRC"] <- replace(m.Screen[ ,t+1 , "CRC"], CRC ==TRUE, 1)
  
  #Assign to CTC, diagnostic colonoscopy or therapeutic colonoscopy depending upon whether adenoma removal/CRC biopsy required
  m.Screen[ ,t+1 , "Attend_CTC"] <- replace(m.Screen[ ,t+1 , "Attend_CTC"], CTC_attend ==TRUE, 1)
  m.Screen[ ,t+1 , "Diagnostic_Scrn_Col"] <- replace(m.Screen[ ,t+1 , "Diagnostic_Scrn_Col"], 
                                                     Col_attend ==TRUE & LR ==FALSE & HR ==FALSE & CRC ==FALSE, 1)
  m.Screen[ ,t+1 , "Therapeutic_Scrn_Col"] <- replace(m.Screen[ ,t+1 , "Therapeutic_Scrn_Col"], LR ==TRUE | HR ==TRUE | CRC ==TRUE, 1)
  
  #Include a repeat diagnostic colonoscopy in inadequates (inadequate for both CTC and colonoscopy)
  m.Screen[ ,t+1 , "Diagnostic_Scrn_Col"] <- m.Screen[ ,t+1 , "Diagnostic_Scrn_Col"] + 
    (((m.Rand[ ,"COL_INAD", t] < inad_COL) & Col_attend) * 1) + (((m.Rand[ ,"CTC_INAD", t] < inad_CTC) & CTC_attend) * 1)
  
  m.Screen
}
