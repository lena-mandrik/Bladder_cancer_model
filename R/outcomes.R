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
  v.cost.diag <- ifelse((m.Diag[, "HG_screen_diag"]==1|m.Diag[, "LG_screen_diag"]==1) & m.Diag[, "yr_diag"]==1, Cost.diag.screen, Cost.diag.sympt)
  v.cost[which(m.Diag[, "yr_diag"]==1)] <- v.cost[which(m.Diag[, "yr_diag"]==1)]+v.cost.diag[which(m.Diag[, "yr_diag"]==1)]
  v.cost
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

aggregate.outcomes <- function(m.M_8s, m.Out, m.M, m.E, m.C, m.Diag, m.Screen, m.State, pop, t) {
  
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
  m.Out["LG_SYMPT", t+1] <- sum((m.Diag[, "LG_new_diag"] ==1 & m.Diag[, "LG_sympt_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St1_SYMPT", t+1] <- sum((m.Diag[, "HG_new_diag"] ==1 & m.Diag[, "HG_sympt_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St2_SYMPT", t+1] <- sum((m.Diag[, "HG_new_diag"] ==1 & m.Diag[, "HG_sympt_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==2) * pop[, "weighting"]) 
  m.Out["HG_St3_SYMPT", t+1] <- sum((m.Diag[, "HG_new_diag"] ==1 & m.Diag[, "HG_sympt_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==3) * pop[, "weighting"]) 
  m.Out["HG_St4_SYMPT", t+1] <- sum((m.Diag[, "HG_new_diag"] ==1 & m.Diag[, "HG_sympt_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==4) * pop[, "weighting"])
  
  m.Out["LG_SCRN", t+1] <- sum((m.Diag[, "LG_new_diag"] ==1 & m.Diag[, "LG_screen_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St1_SCRN", t+1] <- sum((m.Diag[, "HG_new_diag"] ==1 & m.Diag[, "HG_screen_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St2_SCRN", t+1] <- sum((m.Diag[, "HG_new_diag"] ==1 & m.Diag[, "HG_screen_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==2) * pop[, "weighting"]) 
  m.Out["HG_St3_SCRN", t+1] <- sum((m.Diag[, "HG_new_diag"] ==1 & m.Diag[, "HG_screen_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==3) * pop[, "weighting"]) 
  m.Out["HG_St4_SCRN", t+1] <- sum((m.Diag[, "HG_new_diag"] ==1 & m.Diag[, "HG_screen_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==4) * pop[, "weighting"]) 
 
  m.Out["HG_St1_MORT", t+1] <- sum((m.M_8s[, t+1] ==8 & m.M_8s[, t] <8 & m.Diag[, "HG_stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St2_MORT", t+1] <- sum((m.M_8s[, t+1] ==8 & m.M_8s[, t] <8 & m.Diag[, "HG_stage_diag"] ==2) * pop[, "weighting"]) 
  m.Out["HG_St3_MORT", t+1] <- sum((m.M_8s[, t+1] ==8 & m.M_8s[, t] <8 & m.Diag[, "HG_stage_diag"] ==3) * pop[, "weighting"]) 
  m.Out["HG_St4_MORT", t+1] <- sum((m.M_8s[, t+1] ==8 & m.M_8s[, t] <8 & m.Diag[, "HG_stage_diag"] ==4) * pop[, "weighting"])
  
  #weighted counts of harms. Add mortality because of TURBT
  m.Out["Die_TURBT", t+1] <- sum((m.Screen[, t+1, "Die_TURBT"]) * pop[, "weighting"])
  
  #weighted resource use
  m.Out["Invite_DS", t+1] <- sum(m.Screen[, t+1, "Invite_DS"] * pop[, "weighting"]) # total number of dipstick kits sent
  m.Out["Respond_DS", t+1] <- sum(m.Screen[, t+1, "Respond_DS"] * pop[, "weighting"]) # total number of dipstick kits responders
  m.Out["Positive_DS", t+1] <- sum((m.Screen[, t+1, "Positive_DS"]) * pop[, "weighting"]) # total number of positive dipstick 
  m.Out["Invite_FC", t+1] <- sum((m.Screen[, t+1, "Invite_FC"]) * pop[, "weighting"]) # total number of follow up CF that were tested
  m.Out["Diagnostic_FC", t+1] <- sum((m.Screen[, t+1, "Diagnostic_FC"]) * pop[, "weighting"]) # total number of follow up CF that were tested
  m.Out["TURBT_FS", t+1] <- sum((m.Screen[, t+1, "TURBT_FS"]) * pop[, "weighting"]) # total number of follow up CF that were tested
  
  m.Out["FP", t+1] <- sum(m.Screen[, t+1, "FP"] * pop[, "weighting"]) # total number of CTCs used
  m.Out["FN", t+1] <- sum(m.Screen[, t+1, "FN"] * pop[, "weighting"]) # total number of CTCs used
  
  #calculate results per person at model start
  m.Out[, t+1] <- m.Out[, t+1] / sum(pop[, "weighting"])
  
  #Apply discounting to costs, QALYs and Life Years
  costs <- c("TOTAL_COSTS", "BC_COSTS", "SCREEN_COSTS")
  m.Out[costs, t+1] <- m.Out[costs, t+1] * c(rep(v.dwc[t+1], length(costs)))
  m.Out["QALYS", t+1] <- m.Out["QALYS", t+1] * v.dwe[t+1]
  m.Out["LYS", t+1] <- m.Out["LYS", t+1] * v.dwe[t+1]
  
  m.Out
}