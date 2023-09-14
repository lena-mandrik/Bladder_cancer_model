# Define a function to initialize and set up matrices
#' @details
#' This function initialize and set up matrices in the run_simulation script
#' @params
#' nsample: size of the simulated population
#' n.t: number of cycles to run
#' n.s_long: number of states in the matrix with different cancer stages
#' states_long: name of the states
#' run_mode: run mode of the model

initialize_matrices <- function(nsample, n.t, n.s_long, states_long, run_mode) {

  m.M <- m.M_8s <- m.M_8s_KC <- m.C <- m.E <-  matrix(nrow = nsample, ncol = n.t + 1)
  
  # initial state is 1 if 30-yo cohort is run and actual state saved in pop matrix if HSE population is used
  if(cohort==1){
    m.M[, 1] <- m.M_8s[, 1] <- rep(1, nsample)   # indicate the initial health state
  } else if (cohort==0) {
    m.M_8s[, 1] <- m.M[, 1] <- pop[ ,"state"] #if the current pop in England is used, assign the first health state as the current state
    m.M[m.M[, 1]>2 & m.M[, 1] !=4, 1] <- 3 # replace in the short matrix the states for undiagnosed HG cancer as state 3 (all HG BC). The diagnosed HG cancer will be gathered in "DeathOC" state
    
  }
  
  if(run_mode != "Calibration" ){
    m.C[, 1] <- 0 # estimate costs per individual for the initial health state 
    m.E[, 1] <- pop[, "EQ5D"] - (Utility.age * (pop[, "age"] - pop[, "age_0"])) # estimate QALY per individual for the initial health state and starting age
    m.E[, 1] <- replace(m.E[, 1], m.E[, 1] >1, 1) # if EQ-5D over 1 reset to 1
    m.E[, 1] <- replace(m.E[, 1], m.E[, 1] <=(-0.594), -0.594) # if EQ-5D under -0.594 reset to -0.594
    
  }
  
  # create another matrix capturing the current health state for each individual
  m.State <- m.State.KC <-matrix(0, nrow = nsample, ncol = n.s_long)
  colnames(m.State) <- states_long
  for(n in 1:n.s_long) {
    m.State[, n] <- replace(m.State[, n], m.M_8s[, 1] ==n, 1)
    m.State.KC[, n] <- replace(m.State.KC[, n], m.M_8s_KC[, 1] ==n, 1)
    
  }
 
  #Create another matrix for current diagnostic information
  # m.Diag is one as only either BC or KC is assumed in one individual
  m.Diag <- matrix(0, nrow = nsample, ncol = 20)
  
  # When BC said, it means HG
  
  #colnames(m.Diag) <- c("HG_state", "LG_state", "KC_state", "HG_diag", "LG_diag", "HG_age_diag", "LG_age_diag", "HG_sympt_diag", 
  #                      "LG_sympt_diag", "HG_screen_diag", "LG_screen_diag", "FP", 
  # "HG_stage_diag", "yr_diag", "HG_yr_diag", "LG_yr_diag", "HG_yr_onset", "HG_age_onset", "LG_age_onset", "age_BC_death")
  
  
  colnames(m.Diag) <- c("HG_state", "LG_state", "KC_state", "Diag", "LG_diag", "Age_diag", "LG_age_diag", "Sympt_diag", 
                        "LG_sympt_diag", "Screen_diag", "LG_screen_diag", "FP_BC", "FP_KC", 
                        "Stage_diag", "yr_diag", "LG_yr_diag", "yr_onset", "age_onset", "LG_age_onset", "age_C_death")
  
  
  if(run_mode != "Calibration" ){
    # Create an array to gather screening and surveillance information for each cycle
    m.Screen <- array(data = 0, dim = c(nsample, n.t + 1, length(screen_names)))
    dimnames(m.Screen)[[3]] <- screen_names # note that this is defined in set_params
    
    #Create a matrix for gathering aggregate outcomes
    m.Out <- matrix(0, nrow = length(out_names), ncol = n.t + 1)
    rownames(m.Out) <- out_names
    colnames(m.Out) <- c(0:n.t)
  }
  
  #Create additional matrices for subgroup results
  #m.Out_M <- m.Out_F <- m.Out_smoke <- m.Out_past.smoke <- m.Out
  result <- list(
    m.M = m.M,
    m.M_8s = m.M_8s,
    m.M_8s_KC =m.M_8s_KC,
    m.C =m.C,
    m.E =m.E,
    m.Diag=m.Diag,
    m.State =m.State,
    m.Screen =m.Screen,
    m.Out=m.Out
  )
  
  return(result)
}




### Functions to calculate outcomes such as costs and QALYs, which depend upon current health states and age ###

#' @details
#' This function calculates BC treatment costs
#' @params
#' m.State: a matrix of current health states for all individuals
#' m.Diag: a matrix giving diagnostic status for all individuals
#' pop: a population matrix of individuals, each with a current age
#' m.Cost.treat: a matrix of treatment costs by age and stage

#' @return a vector of costs for each person

f.calc.cost <- function(m.State, m.Diag, m.Cost.treat){
  
  #Treatment costs in people diagnosed with cancer (up to 5 years after diagnosis) 
  m.cost <- (m.State * m.Cost.treat[paste(m.Diag[, "yr_diag"]), ]) 
  
  #Add costs of diagnosis to symptomatically diagnosed only (screen diagnosed are added on a later stage in aggregated outcomes)
  v.cost.diag <- ifelse((m.Diag[, "HG_sympt_diag"]==1|m.Diag[, "LG_sympt_diag"]==1) & m.Diag[, "yr_diag"]==1, Cost.diag.sympt, 0)
  
  # bind the diagnostic costs to treatment costs
  m.cost <- cbind(m.cost, v.cost.diag)
    
  v.cost <-rowSums(m.cost)
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
f.calc.utility <- function(m.State, m.Diag, pop, t) {    
  
  #First calculate age decrement
  v.utility <- pop[, "EQ5D"] - (Utility.age * (pop[, "age"] - pop[, "age_0"]))
  
  # Mark those individuals with the diagnosed disease
  LG_u <- (m.Diag[, "LG_yr_diag"] >0 & m.Diag[, "LG_yr_diag"] <=3 & m.State[,"BC_LG"] ==1)*1 # Impact on utilities for LG cancers is during 3 years
  HG_1.3_u <- (m.Diag[, "HG_yr_diag"] >0 & m.Diag[, "HG_yr_diag"] <=5 & (m.State[,"St1_HG"] ==1 |m.State[,"St2_HG"] ==1 | m.State[,"St3_HG"] ==1))*1 # Impact on utilities for HG cancers is during 5 years
  HG_4_u <- (m.Diag[, "HG_yr_diag"] >0 & m.Diag[, "HG_yr_diag"] <=5 & m.State[,"St4_HG"] ==1)*1 # Impact on utilities for HG cancers is during 5 years
  
  # Assign disutility for each diagnosed state by using the multipliers
  disutility_BC <- LG_u*Disutility.LG + HG_1.3_u*Disutility.HG.St1.3 + HG_4_u*Disutility.HG.St4
  
  #substract cancer disutility from the individual EQ5D in people with EQ5D >= -0.594. Assume the utilities in people with lower utilities remain so
  v.utility[v.utility >=-0.594] <- (v.utility*disutility_BC)[v.utility >=-0.594]
  
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

f.aggregate.outcomes <- function(m.M_8s, m.Out, m.M, m.E, m.C, m.Diag, m.Screen, m.State, pop, t) {
  
  ###Fill outcomes matrix
  
  #weighted costs of BC treatment, half cycle corrected
  m.Out["BC_COSTS", t+1] <- sum(0.5 * (m.C[, t] + m.C[, t+1]) * pop[, "weighting"]) # total weighted CRC treatment costs
  
  #weighted costs of DT screening (excluding follow-up), half cycle corrected
  m.Out["SCREEN_COSTS", t+1] <- sum(0.5 * ((m.Screen[, t+1, DS_names] %*% m.Cost.screen[DS_names, ]) + (m.Screen[, t, DS_names] %*% m.Cost.screen[DS_names, ])) * pop[, "weighting"]) 
  
  m.Out["DIAG_COSTS", t+1] <- sum(0.5*((m.Screen[, t+1, "Respond_diag"] * Cost.diag.screen1) + (m.Screen[, t+1, "Respond_Cyst"] * Cost.diag.screen2))+
                                  0.5*((m.Screen[, t, "Respond_diag"] * Cost.diag.screen1) + (m.Screen[, t, "Respond_Cyst"] * Cost.diag.screen2)))
  
  #total weighted costs, half cycle corrected
  m.Out["TOTAL_COSTS", t+1] <- m.Out["BC_COSTS", t+1] + m.Out["SCREEN_COSTS", t+1] + m.Out["DIAG_COSTS", t+1]
  
  #weighted QALYs and life years, half cycle corrected
  m.Out["QALYS", t+1] <- sum(0.5 * (m.E[, t] + m.E[, t+1]) * pop[, "weighting"]) # total weighted QALYs
  m.Out["LYS", t+1] <- sum(0.5 * ((m.M[, t] < 4) + (m.M[, t+1] < 4)) * pop[, "weighting"]) # total weighted LYs
  
  #weighted cancer diagnoses and mortality
  m.Out["LG_SYMPT", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "LG_sympt_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St1_SYMPT", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "HG_sympt_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St2_SYMPT", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "HG_sympt_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==2) * pop[, "weighting"]) 
  m.Out["HG_St3_SYMPT", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "HG_sympt_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==3) * pop[, "weighting"]) 
  m.Out["HG_St4_SYMPT", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "HG_sympt_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==4) * pop[, "weighting"])
  
  m.Out["LG_SCRN", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "LG_screen_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St1_SCRN", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "HG_screen_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==1) * pop[, "weighting"]) 
  m.Out["HG_St2_SCRN", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "HG_screen_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==2) * pop[, "weighting"]) 
  m.Out["HG_St3_SCRN", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "HG_screen_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==3) * pop[, "weighting"]) 
  m.Out["HG_St4_SCRN", t+1] <- sum((m.Diag[ ,"yr_diag"]==1 & m.Diag[, "HG_screen_diag"] ==1 & m.Diag[, "HG_stage_diag"] ==4) * pop[, "weighting"]) 
 
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
  m.Out["Respond_Cyst", t+1] <- sum((m.Screen[, t+1, "Respond_Cyst"]) * pop[, "weighting"]) # total number of follow up CF that were tested
  m.Out["Diagnostic_Cyst", t+1] <- sum((m.Screen[, t+1, "Diagnostic_Cyst"]) * pop[, "weighting"]) # total number of follow up CF that were tested
  m.Out["TURBT", t+1] <- sum((m.Screen[, t+1, "TURBT"]) * pop[, "weighting"]) # total number of follow up CF that were tested
  
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