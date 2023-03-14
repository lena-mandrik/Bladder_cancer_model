### Functions to calculate outcomes such as costs and QALYs, which depend upon current health states and age ###

#' @details
#' This function calculates CRC treatment costs
#' @params
#' m.State: a matrix of current health states for all individuals
#' m.Diag: a matrix giving diagnostic status for all individuals
#' pop: a population matrix of individuals, each with a current age
#' m.Cost.treat: a matrix of treatment costs by age and stage

#' @return a vector of costs for each person

calc.cost <- function(m.State, m.Diag, m.Cost.treat){
  
  #Treatment costs in people diagnosed with cancer (up to 5 years after diagnosis) 
  m.cost <- (m.State * m.Cost.treat[paste(m.Diag[, "yr_diag"]), ]) 
  v.cost <-rowSums(m.cost)
  
  #Add costs of diagnosis depending on whether screen- or symptomatic diagnosed
  v.cost.diag <- ifelse(m.Diag[, "Screen_diag"]==1 & m.Diag[, "yr_diag"]==1, Cost.diag.screen, Cost.diag.sympt)
  v.cost[which(m.Diag[, "yr_diag"]==1)] <- v.cost[which(m.Diag[, "yr_diag"]==1)]+v.cost.diag
}

# !!! currently the impact on costs of smoking and no smoking is not implemented

############################################################################################################

#' @details
#' This function calculates utilities
#' @params
#' m.State: a matrix of current health states for all individuals
#' m.Diag: a matrix giving diagnostic status for all individuals
#' m.Screen: an array of screening history
#' pop: a population matrix of individuals, each with a current age
#' t: cycle
#' @return a vector of utilities for each person
calc.utility <- function(m.State, m.Diag, pop, t) {    
  
  #First calculate age decrement
  v.utility <- pop[, "EQ5D"] - (Utility.age * (pop[, "age"] - pop[, "age_0"]))
  
  # Mark those individuals with the diagnosed disease
  LG_u <- (m.Diag[, "LG_yr_diag"] >0 & m.Diag[, "LG_yr_diag"] <3)*1 # Impact on utilities for LG cancers is during 3 years
  HG_1.3_u <- (m.Diag[, "HG_yr_diag"] >0 & m.Diag[, "HG_yr_diag"] <5 & (m.State[,"St1_HG"] ==1 |m.State[,"St2_HG"] ==1 | m.State[,"St3_HG"] ==1))*1 # Impact on utilities for HG cancers is during 5 years
  HG_4_u <- (m.Diag[, "HG_yr_diag"] >0 & m.Diag[, "HG_yr_diag"] <5 & m.State[,"St4_HG"] ==1)*1 # Impact on utilities for HG cancers is during 5 years
  
  # Assign disutility for each diagnosed state
  disutility_BC <- LG_u*Disutility.LG + HG_1.3_u*Disutility.HG.St1.3 + HG_4_u*Disutility.HG.St4
  
  #substract cancer disutility from the individual EQ5D in people with EQ5D >= -0.594. Assume the utilities in people with lower utilities remain so
  v.utility[v.utility >=-0.594] <- (v.utility-disutility_BC)[v.utility >=-0.594]
  
  #Subtract any transient utility decrements associated with screening (note harm_names defined in set_params)
  #v.utility <- v.utility + (m.Screen[, t+1, ] %*% m.Utility_harms)
  
  #Replace any individual EQ5D score below minimum value, with minimum (-0.594)
  v.utility <- replace(v.utility, v.utility < -0.594, -0.594)
  
  #Replace any individual EQ5D score above maximum value, with maximum (1)
  v.utility <- replace(v.utility, v.utility > 1, 1)
  
  #Replace EQ5D score with 0 for dead people
  v.utility <- replace(v.utility, m.State[, "DeathOC"] ==1 | m.State[, "DeathBC"] ==1, 0)
  
  v.utility
}  


#' @details
#' This function calculates aggregate, weighted outcomes for each year of simulation, with half cycle correction for QALYs and Life years
#' @params
#' m.Out: an outcomes matrix that is empty for time t.
#' m.M: a matrix of health states for all individuals over time
#' m.C: a matrix of costs for all individuals over time
#' m.E: a matrix of utilities for all individuals over time
#' m.Diag: a matrix giving diagnostic status for all individuals
#' m.Screen: an array of screening history
#' m.State: a matrix of current health states for all individuals
#' pop: a population matrix of individuals, each with a current age
#' t: cycle number
#' @return a matrix of aggregated, weighted outcomes for time t

aggregate.outcomes <- function(m.Out, m.M, m.E, m.C, m.Diag, m.Screen, m.State, pop, t) {
  
  ###Fill outcomes matrix
  
  #weighted costs of CRC treatment, half cycle corrected
  m.Out["BC_COSTS", t+1] <- sum(0.5 * (m.C[, t] + m.C[, t+1]) * pop[, "weighting"]) # total weighted CRC treatment costs
  
  #weighted costs of DT screening (excluding follow-up), half cycle corrected
  m.Out["SCREEN_COSTS", t+1] <- sum(0.5 * ((m.Screen[, t+1, DS_names] %*% m.Cost.screen[DS_names, ]) + (m.Screen[, t, DS_names] %*% m.Cost.screen[DS_names, ])) * pop[, "weighting"]) 
  
  #total weighted costs, half cycle corrected
  m.Out["TOTAL_COSTS", t+1] <- m.Out["BC_COSTS", t+1] + m.Out["SCREEN_COSTS", t+1] 
  
  #weighted QALYs and life years, half cycle corrected
  m.Out["QALYS", t+1] <- sum(0.5 * (m.E[, t] + m.E[, t+1]) * pop[, "weighting"]) # total weighted QALYs
  m.Out["LYS", t+1] <- sum(0.5 * ((m.M[, t] < 4) + (m.M[, t+1] < 4)) * pop[, "weighting"]) # total weighted LYs
  
  #weighted cancer diagnoses and mortality
  m.Out["LG_SYMPT", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "sympt_diag"] ==1 & m.Diag[, "stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St1_SYMPT", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "sympt_diag"] ==1 & m.Diag[, "stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St2_SYMPT", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "sympt_diag"] ==1 & m.Diag[, "stage_diag"] ==2) * pop[, "weighting"]) 
  m.Out["HG_St3_SYMPT", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "sympt_diag"] ==1 & m.Diag[, "stage_diag"] ==3) * pop[, "weighting"]) 
  m.Out["HG_St4_SYMPT", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "sympt_diag"] ==1 & m.Diag[, "stage_diag"] ==4) * pop[, "weighting"]) 
  m.Out["LG_SCRN", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "surv_diag"] ==1 & m.Diag[, "stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St1_SCRN", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "screen_diag"] ==1 & m.Diag[, "stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St2_SCRN", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "screen_diag"] ==1 & m.Diag[, "stage_diag"] ==2) * pop[, "weighting"]) 
  m.Out["HG_St3_SCRN", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "screen_diag"] ==1 & m.Diag[, "stage_diag"] ==3) * pop[, "weighting"]) 
  m.Out["HG_St4_SCRN", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "screen_diag"] ==1 & m.Diag[, "stage_diag"] ==4) * pop[, "weighting"]) 
  #m.Out["CRC_B_SURV", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "surv_diag"] ==1 & m.Diag[, "stage_diag"] ==2) * pop[, "weighting"]) 
  #m.Out["CRC_C_SURV", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "surv_diag"] ==1 & m.Diag[, "stage_diag"] ==3) * pop[, "weighting"]) 
  #m.Out["CRC_D_SURV", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "surv_diag"] ==1 & m.Diag[, "stage_diag"] ==4) * pop[, "weighting"]) 
  #m.Out["CRC_A_INT", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "interval_diag"] ==1 & m.Diag[, "stage_diag"] ==1) * pop[, "weighting"]) 
  #m.Out["CRC_B_INT", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "interval_diag"] ==1 & m.Diag[, "stage_diag"] ==2) * pop[, "weighting"]) 
  #m.Out["CRC_C_INT", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "interval_diag"] ==1 & m.Diag[, "stage_diag"] ==3) * pop[, "weighting"]) 
  #m.Out["CRC_D_INT", t+1] <- sum((m.Diag[, "new_diag"] ==1 & m.Diag[, "interval_diag"] ==1 & m.Diag[, "stage_diag"] ==4) * pop[, "weighting"]) 
  m.Out["CRC_A_MORT", t+1] <- sum((m.M[, t+1] ==8 & m.M[, t] <8 & m.Diag[, "stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["CRC_B_MORT", t+1] <- sum((m.M[, t+1] ==8 & m.M[, t] <8 & m.Diag[, "stage_diag"] ==2) * pop[, "weighting"]) 
  m.Out["CRC_C_MORT", t+1] <- sum((m.M[, t+1] ==8 & m.M[, t] <8 & m.Diag[, "stage_diag"] ==3) * pop[, "weighting"]) 
  m.Out["CRC_D_MORT", t+1] <- sum((m.M[, t+1] ==8 & m.M[, t] <8 & m.Diag[, "stage_diag"] ==4) * pop[, "weighting"])
  
  #weighted counts of harms
  m.Out["HARM_BLEED", t+1] <- sum((m.Screen[, t+1, "Bleed_FS"] + m.Screen[, t+1, "Bleed_Colon"]) * pop[, "weighting"]) 
  m.Out["HARM_PERF", t+1] <- sum((m.Screen[, t+1, "Perf_FS"] + m.Screen[, t+1, "Perf_Colon"] + m.Screen[, t+1, "Perf_CTC"]) * pop[, "weighting"])
  m.Out["HARM_MORT", t+1] <- sum((m.Screen[, t+1, "Die_FS"] + m.Screen[, t+1, "Die_Colon"] + m.Screen[, t+1, "Die_CTC"]) * pop[, "weighting"])
  
  #weighted resource use
  m.Out["FS_USE", t+1] <- sum((m.Screen[, t+1, "Diagnostic_FS"] + m.Screen[, t+1, "Therapeutic_FS"]) * pop[, "weighting"]) # total number of FS usage
  m.Out["FIT_INVITE", t+1] <- sum(m.Screen[, t+1, "Invite_FIT"] * pop[, "weighting"]) # total number of FIT invites sent
  m.Out["FIT_RESPOND", t+1] <- sum(m.Screen[, t+1, "Respond_FIT"] * pop[, "weighting"]) # total number of FIT responders
  m.Out["SCREEN_COL_USE", t+1] <- sum((m.Screen[, t+1, "Diagnostic_Scrn_Col"] + m.Screen[, t+1, "Therapeutic_Scrn_Col"]) * pop[, "weighting"]) # total number of screening colonoscopies used
  m.Out["SURV_COL_USE", t+1] <- sum((m.Screen[, t+1, "Diagnostic_Surv_Col"] + m.Screen[, t+1, "Therapeutic_Surv_Col"]) * pop[, "weighting"]) # total number of surveillance colonoscopies used
  m.Out["CTC_USE", t+1] <- sum(m.Screen[, t+1, "Attend_CTC"] * pop[, "weighting"]) # total number of CTCs used
  
  #calculate results per person at model start
  m.Out[, t+1] <- m.Out[, t+1] / sum(pop[, "weighting"])
  
  #Apply discounting to costs, QALYs and Life Years
  costs <- c("TOTAL_COSTS", "CRC_COSTS", "FIT_COSTS", "FS_COSTS", "CTC_COSTS", "SCREEN_COL_COSTS", "SURV_COSTS", "HARM_COSTS")
  m.Out[costs, t+1] <- m.Out[costs, t+1] * c(rep(v.dwc[t+1], length(costs)))
  m.Out["QALYS", t+1] <- m.Out["QALYS", t+1] * v.dwe[t+1]
  m.Out["LYS", t+1] <- m.Out["LYS", t+1] * v.dwe[t+1]
  
  m.Out
}