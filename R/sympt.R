
#################################################################################################################################
# Function assigning people to have recorded symptomatic diagnosis if they have cancer based on probabilities
#################################################################################################################################
#' @details
#' This function updates m.Diag following symptomatic detection of CRC
#' @params
#' m.Diag: a matrix giving diagnostic status for all individuals
#' m.State: a matrix of current health states for all individuals
#' m.Rand: a random number array
#' pop: a population matrix of individuals, each with a current age
#' t: current time point
#' @return an updated current diagnostic information matrix


f.symptom <- function(m.Diag, m.State, m.Rand, pop, t, m.M, elig_time, nsample) {
  

  #Specify who is eligible for diagnosis of those still alive (those not yet diagnosed)
  elig_HG <- m.Diag[, "Diag"] ==0 & m.M[ ,t+1] ==3 &  elig_time ==1 # for HG cancer
  elig_LG <- m.Diag[, "LG_diag"]==0 & m.M[ ,t+1]==2 &  elig_time ==1 # for LG cancer
  elig_KC <- m.Diag[, "Diag"] ==0 & m.M[ ,t+1] ==4 &  elig_time ==1 # for HG cancer
  
  # Probability to become symptomatic patient
  # Returns an annual probability to be diagnosed by time of cancer onset if a person has cancer
  
  
  #Update m.Diag matrix for symptomatic diagnosis for high grade cancers
  #Determine who is newly diagnosed
  New_HG <- 1*((m.Rand[ ,"SYMPT_HG", t] < risk.Sympt.diag_HG) & elig_HG)
  m.Diag[, "HG_yr_diag"] <- m.Diag[, "HG_yr_diag"] + New_HG #m.Diag[, "HG_new_diag"]
  m.Diag[, "HG_diag"] <- m.Diag[, "HG_diag"] + New_HG #m.Diag[, "HG_new_diag"]
  m.Diag[, "HG_sympt_diag"] <- m.Diag[, "HG_sympt_diag"] + New_HG #m.Diag[, "HG_new_diag"]
  m.Diag[, "HG_age_diag"] <- m.Diag[, "HG_age_diag"] + (pop[, "age"] * New_HG) #m.Diag[, "HG_new_diag"])
  m.Diag[, "HG_stage_diag"] <- m.Diag[, "HG_stage_diag"] + 
    ((m.State[, "St1_HG"] * 1 + m.State[, "St2_HG"] * 2 + m.State[, "St3_HG"] * 3+ m.State[, "St4_HG"] * 4) * New_HG)
  
  #Update m.Diag matrix for symptomatic diagnosis for low grade cancers
  New_LG <- 1*((m.Rand[ ,"SYMPT_LG", t] < risk.Sympt.diag_LG) & elig_LG)
  m.Diag[, "LG_yr_diag"] <- m.Diag[, "LG_yr_diag"] + New_LG#m.Diag[, "LG_new_diag"]
  m.Diag[ , "LG_diag"] <- m.Diag[ , "LG_diag"] + New_LG#m.Diag[, "LG_new_diag"]
  m.Diag[ , "LG_age_diag"] <- m.Diag[, "LG_age_diag"] + (pop[, "age"] * New_LG)#m.Diag[, "LG_new_diag"])
  m.Diag[, "LG_sympt_diag"] <- m.Diag[, "LG_sympt_diag"] + New_LG #m.Diag[, "LG_new_diag"]

  m.Diag[, "yr_diag"] <- m.Diag[, "LG_yr_diag"]
  m.Diag[, "yr_diag"] <- replace(m.Diag[, "yr_diag"], m.Diag[, "HG_diag"] >0, m.Diag[m.Diag[, "HG_diag"] >0, "HG_yr_diag"])
  
  m.Diag
  
} # close the function bracket 

#################################################################################################################################
# Function calculating the risk of symptomatic diagnosis for either BCHG, BCLG, or KC
#################################################################################################################################
#' @details
#' This function updates m.Diag following symptomatic detection of CRC
#' @params
#' P.sympt.diag_Age: a coefficient incorporating a decrease in symptomatic presentation by age
#' risk.Sympt.diag: a matrix with calibrated parameters on symptomatic recording rate
#' pop: a population matrix of individuals, each with a current age
#' 
#' @return an updated risk.Sympt.diag

# Update for those who are older than 80 with the age parameter for a symptomatic presentation rate for HGBC stage 1,2 and LGBC

f.risk.Sympt.diag <- function(nsample, Cancer.env, m.State){
 
  age_sympt <- cbind(pop[, "age"],rep(1, nsample))
  
  age_sympt[age_sympt[,1]>=75, 2] <- age_sympt[age_sympt[,1]>=75, 2] * (Cancer.env$P.sympt.diag_Age ^ (age_sympt[age_sympt[,1]>=75, 1] - 75))
  
  # risk to become symptomatic is calibrated for each stage
  risk.Sympt.diag <- (m.State %*% Cancer.env$Symp_params)*age_sympt[,2]
  risk.Sympt.diag[risk.Sympt.diag>=1] <- 0.999 #replace the individual transitions of diagnosis if exceed 1 with 0.99 probability to be diagnosed
  
 return(risk.Sympt.diag)
  
}

