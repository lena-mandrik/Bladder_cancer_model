#' @details
#' This function updates m.Diag following symptomatic detection of CRC
#' @params
#' m.Diag: a matrix giving diagnostic status for all individuals
#' m.State: a matrix of current health states for all individuals
#' m.Rand: a random number array
#' pop: a population matrix of individuals, each with a current age
#' t: current time point
#' @return an updated current diagnostic information matrix

f.symptom <- function(m.Diag, m.State, m.Rand, pop, t, m.M) {
  
  #Reset newly diagnosed and move year of diagnosis on by one
  m.Diag[, "New_diag"] <- 0
  new_diag_LG <- 0
  #m.Diag[, "HG_yr_diag"][m.Diag[, "HG_yr_diag"] >=1] <- m.Diag[, "HG_yr_diag"][m.Diag[, "HG_yr_diag"] >=1] + 1
  #m.Diag[, "LG_yr_diag"][m.Diag[, "LG_yr_diag"] >=1] <- m.Diag[, "LG_yr_diag"][m.Diag[, "LG_yr_diag"] >=1] + 1
  m.Diag[, "yr_diag"][m.Diag[, "yr_diag"] >=1] <- m.Diag[, "yr_diag"][m.Diag[, "yr_diag"] >=1] + 1
  
  #Specify who is eligible for diagnosis of those still alive (those not yet diagnosed)
  elig_HG <- (m.Diag[, "BC_diag"] ==0 & m.M[ ,t+1] ==3)||(m.M[ ,t+1] ==3 & m.Diag[, "HG_stage_diag"]==0 & m.Diag[, "LG_state"]==1) # for HG cancer either those who hasn't been diagnosed or those who had been diagnosed with LG and progressed to HG
  elig_LG <- m.Diag[, "BC_diag"]==0 & m.M[ ,t+1]==2 # for LG cancer
  
  # Probability to become symptomatic patient
  # m.Diag[ ,"BC_state"] - should be 1 for each patient with BC in any HG state.
  # Returns an annual probability to be diagnosed by time of cancer onset if a person has cancer
  
  # TP to become symptomatic (TP.Sympt.diag) is calculated from two calibrated parameters - P.sympt.diag_A_HGBC and P.sympt.diag_B_HGBC.
  # This is because the symptomatic presentation rate is not linear and increase with the time since onset similar to the stage by TNM
  
  TP.Sympt.diag <- m.Diag[ ,"HG_state"]*P.sympt.diag_A_HGBC*P.sympt.diag_B_HGBC^m.Diag[ ,"HG_yr_onset"]
  
  # Update for those who are older than 80 with the decrement in a symptomatic presentation rate for HGBC
  age_sympt_decrement <- P.sympt.diag_Age80_HGBC^(pop[, "age"]- 80)
  age_sympt_decrement <- replace(age_sympt_decrement, age_sympt_decrement>1,1)
  risk.Sympt.diag_HG <- TP.Sympt.diag * age_sympt_decrement
  
  # Update for those who are older than 80 with the decrement in a symptomatic presentation rate for LGBC
  
  #age_sympt_decrement_LG <-P.sympt.diag_LGBC^(pop[, "age"]- 80)
  #age_sympt_decrement_LG <- replace(age_sympt_decrement_LG, age_sympt_decrement_LG>1,1)
  risk.Sympt.diag_LG <- age_sympt_decrement * P.sympt.diag_LGBC
  
  #Determine who is newly diagnosed
  new_diag_HG <- (m.Rand[ ,"SYMPT_HG", t] < risk.Sympt.diag_HG) & elig_HG #for HG cancer
  new_diag_LG <- (m.Rand[ ,"SYMPT_LG", t] < risk.Sympt.diag_LG) & elig_LG #for LG cancer
  
  #Update m.Diag matrix for symptomatic diagnosis for high grade cancers
  m.Diag[, "New_diag"] <- new_diag_HG + new_diag_LG
  m.Diag[, "yr_diag"] <- m.Diag[, "yr_diag"] + m.Diag[, "New_diag"]
  #m.Diag[, "HG_diag"] <- m.Diag[, "HG_diag"] + m.Diag[, "HG_new_diag"]
  m.Diag[, "BC_diag"] <- m.Diag[, "BC_diag"] + m.Diag[, "New_diag"]
  #m.Diag[, "HG_sympt_diag"] <- m.Diag[, "HG_sympt_diag"] + m.Diag[, "HG_new_diag"]
  m.Diag[, "Sympt_diag"] <- m.Diag[, "Sympt_diag"] + m.Diag[, "New_diag"]
  
  m.Diag[, "HG_age_diag"] <- m.Diag[, "HG_age_diag"] + (pop[, "age"] * new_diag_HG)
  #m.Diag[, "HG_age_diag"] <- m.Diag[, "HG_age_diag"] + (pop[, "age"] * m.Diag[, "HG_new_diag"])
  
  m.Diag[, "HG_stage_diag"] <- m.Diag[, "HG_stage_diag"] + 
    ((m.State[, "St1_HG"] * 1 + m.State[, "St2_HG"] * 2 + m.State[, "St3_HG"] * 3+ m.State[, "St4_HG"] * 4) * m.Diag[, "New_diag"]*m.Diag[, "HG_state"])
  
  #Update m.Diag matrix for symptomatic diagnosis for low grade cancers
  #m.Diag[ , "LG_diag"] <- m.Diag[ , "LG_diag"] + new_diag_LG
  m.Diag[ , "LG_age_diag"] <- m.Diag[, "LG_age_diag"] + (pop[, "age"] * new_diag_LG)
  #m.Diag[, "LG_yr_diag"] <- m.Diag[, "LG_yr_diag"] + new_diag_LG
  
  #Update yr_diag so that it is yr since diagnosis for LG cancer or HG if LG progressed to HG
 # m.Diag[, "yr_diag"][m.Diag[, "LG_yr_diag"] >=1] <- m.Diag[, "LG_yr_diag"][m.Diag[, "LG_yr_diag"] >=1]
 # m.Diag[, "yr_diag"][m.Diag[, "HG_yr_diag"] >=1] <- m.Diag[, "HG_yr_diag"][m.Diag[, "HG_yr_diag"] >=1]
  
  m.Diag
  
} # close the function bracket 
