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
  m.Diag[, "HG_new_diag"] <- m.Diag[, "LG_new_diag"] <-0
  m.Diag[, "HG_yr_diag"][m.Diag[, "HG_yr_diag"] >=1] <- m.Diag[, "HG_yr_diag"][m.Diag[, "HG_yr_diag"] >=1] + 1
  m.Diag[, "LG_yr_diag"][m.Diag[, "LG_yr_diag"] >=1] <- m.Diag[, "LG_yr_diag"][m.Diag[, "LG_yr_diag"] >=1] + 1
  
  #Specify who is eligible for diagnosis of those still alive (those not yet diagnosed)
  elig_HG <- m.Diag[, "HG_diag"] ==0 & m.M[ ,t+1] ==3 # for HG cancer
  elig_LG <- m.Diag[, "LG_diag"]==0 & m.M[ ,t+1]==2 # for LG cancer
  
  # Probability to become symptomatic patient
  # Returns an annual probability to be diagnosed by time of cancer onset if a person has cancer
  
  # Update for those who are older than 70 with the age parameter for a symptomatic presentation rate for HGBC
  age_sympt <- cbind(pop[, "age"],rep(1, n.i), rep(1, n.i))
  age_sympt[,2][age_sympt[,1]>=67] <- age_sympt[,2][age_sympt[,1]>=67]*(P.sympt.diag_Age_HGBC^(age_sympt[,1][age_sympt[,1]>=67]- 67))
  
  # risk to become symptomaticis calculated from two calibrated parameters - P.sympt.diag_HGBC and P.sympt.diag_t_HGBC.
  # P.sympt.diag_HGBC is the probability to become symptomatic the 1st year after the onset and P.sympt.diag_t_HGBC is the coefficient increasing the p of symptoms with time since onset
  
  risk.Sympt.diag <- m.Diag[ ,"HG_state"]*P.sympt.diag_HGBC*age_sympt[,2]+(P.sympt.diag_t_HGBC*m.Diag[ ,"HG_yr_onset"])
  
  risk.Sympt.diag[risk.Sympt.diag>=1] <- 0.999 #replace the individual transitions of diagnosis if exceed 1 with 0.99 probability to be diagnosed
  
  # Update for those who are older than 80 with the decrement in a symptomatic presentation rate for LGBC
  
  age_sympt[,3][age_sympt[,1]>=67] <- age_sympt[,3][age_sympt[,1]>=67]*(P.sympt.diag_Age_LGBC^(age_sympt[,1][age_sympt[,1]>=67]- 67))
  
  risk.Sympt.diag_LG <- P.sympt.diag_LGBC*age_sympt[,3]
  risk.Sympt.diag_LG[risk.Sympt.diag_LG>=1] <- 0.999 #replace the individual transitions of diagnosis if exceed 1 with 0.99 probability to be diagnosed
  
  #Update m.Diag matrix for symptomatic diagnosis for high grade cancers
  #Determine who is newly diagnosed
  m.Diag[, "HG_new_diag"] <- 1*((m.Rand[ ,"SYMPT_HG", t] < risk.Sympt.diag) & elig_HG) #for HG cancer
  m.Diag[, "HG_yr_diag"] <- m.Diag[, "HG_yr_diag"] + m.Diag[, "HG_new_diag"]
  m.Diag[, "HG_diag"] <- m.Diag[, "HG_diag"] + m.Diag[, "HG_new_diag"]
  m.Diag[, "HG_sympt_diag"] <- m.Diag[, "HG_sympt_diag"] + m.Diag[, "HG_new_diag"]
  m.Diag[, "HG_age_diag"] <- m.Diag[, "HG_age_diag"] + (pop[, "age"] * m.Diag[, "HG_new_diag"])
  m.Diag[, "HG_stage_diag"] <- m.Diag[, "HG_stage_diag"] + 
    ((m.State[, "St1_HG"] * 1 + m.State[, "St2_HG"] * 2 + m.State[, "St3_HG"] * 3+ m.State[, "St4_HG"] * 4) * m.Diag[, "HG_new_diag"])
  
  #Update m.Diag matrix for symptomatic diagnosis for low grade cancers
  m.Diag[, "LG_new_diag"] <- 1*((m.Rand[ ,"SYMPT_LG", t] < risk.Sympt.diag_LG) & elig_LG) #for LG cancer
  m.Diag[, "LG_yr_diag"] <- m.Diag[, "LG_yr_diag"] + m.Diag[, "LG_new_diag"]
  m.Diag[ , "LG_diag"] <- m.Diag[ , "LG_diag"] + m.Diag[, "LG_new_diag"]
  m.Diag[ , "LG_age_diag"] <- m.Diag[, "LG_age_diag"] + (pop[, "age"] * m.Diag[, "LG_new_diag"])
  m.Diag[, "LG_sympt_diag"] <- m.Diag[, "LG_sympt_diag"] + m.Diag[, "LG_new_diag"]
  m.Diag[, "LG_age_diag"] <- m.Diag[, "LG_age_diag"] + (pop[, "age"] * m.Diag[, "LG_new_diag"])
  
  m.Diag[, "yr_diag"] = m.Diag[, "LG_yr_diag"]
  m.Diag[, "yr_diag"] = replace(m.Diag[, "yr_diag"], m.Diag[, "HG_diag"] >0, m.Diag[m.Diag[, "HG_diag"] >0, "HG_yr_diag"])
  
  m.Diag
  
} # close the function bracket 
