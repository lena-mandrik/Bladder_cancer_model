### Functions for screening ###

#' @details
#' This function calculates individualised screening parameters for each person
#' @params
#' pop: a population matrix of individuals, each with a current age
#' m:Screen: an array of screening history
#' m.State: a matrix giving current health state for all individuals
#' @return a matrix of individualised screening parameters
#' 

f.calc.screen.Params <- function(pop, m.Screen, m.State) {
  
    #Calculates dipstick uptake by personal characteristics
    DS_uptake <- 1/(1+exp(-((cbind((m.State[, "DeathBC"] ==0 & m.State[, "DeathOC"] ==0), pop[,"age"] <55, 
                                    pop[,"age"] >=55 & pop[,"age"] <60, pop[,"age"] >=65 & pop[,"age"] <70, pop[,"age"] >=70, 
                                    pop[, "sex"] ==0, (rowSums(m.Screen[, ,"Invite_DS"]) >=1 & rowSums(m.Screen[, ,"Respond_DS"]) ==0),
                                    (rowSums(m.Screen[, ,"Respond_DS"]) >=1), pop[, "imd"] ==2, pop[, "imd"] ==3, pop[, "imd"] ==4, 
                                    pop[, "imd"] ==5, pop[, "ethnic"] ==3) * 1) %*% coef_DT_Uptk)))
    
    #Calculates DT sensitivity or false positives depending upon underlying health state
    DS_diag <- (m.State %*% test_accuracy[, "Sens"]) 
    
    #Add diagnostic uptake: considered as 1 in the basecase analysis
    Diag_uptake <- rep(Diag.UPTK, n.i)
    
    scr.Params <- cbind(DS_uptake, DS_diag, Diag_uptake)
    
    colnames(scr.Params) <- c("DS_uptake", "DS_diag", "Cyst_uptake")
  
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
#' DS_age - age when screening starts
#' DS_round - number of the screening rounds
#' DS_freq - frequency of the screening rounds
#' @return an updated screening history array
#' 
f.DS_screen <- function(m.Screen, m.Diag, m.State, m.Rand, pop, t, scr.Params, DS_age, DS_round, DS_freq, n_round) {
  
  # Determine who is eligible for screening
  #m.Screen[ ,t+1 , "Invite_DS"] <- (m.Diag[, "LG_diag"] ==0 & m.Diag[, "HG_diag"] ==0 & 
  #  m.State[, "DeathBC"] ==0 & m.State[, "DeathOC"] ==0 &
  #  pop[, "age"] == DS_age)*1
  
  # set n_round as the current screening round
  
  # set eligibility by smoking status (considering that smoking status is variable)
  smoke_eligible <- as.matrix(rep(1, n.i), ncol=1) #set eligibility of the population by their current smoking status on that year
  if(screen_elig == "c.smoke"){
    smoke_eligible[pop[ ,"current_smoke"] !=1, ] <- 0
    }else if(screen_elig == "all.smoke"){
      smoke_eligible[pop[ ,"current_smoke"] !=1 | pop[ ,"past_smoke"] !=1 , ] <- 0
    }
  
  m.Screen[ ,t+1 , "Invite_DS"] <- (m.Diag[, "LG_diag"] ==0 & m.Diag[, "HG_diag"] ==0 & # population not diagnosed with cancer yet
                                      m.State[, "DeathBC"] ==0 & m.State[, "DeathOC"] ==0 & # alive population
                                      smoke_eligible[,1]==1 & #the pop status by their current smoking eligibility
                                      (pop[, "age"] == DS_age || # for one-time screening, should be equal to the pop age; or the screening start for the repetitive screening
                                         (pop[, "age"] > DS_age & n_round < DS_round & (DS_freq ==1 | # for annual screening
                                            (m.Screen[ ,t, "Invite_DS"]==0 & DS_freq ==2)))))*1 # for biennial screening
  
  # Update the number of screen rounds each person received 
  n_round[m.Screen[ ,t+1 , "Invite_DS"]==1] <<- n_round[m.Screen[ ,t+1 , "Invite_DS"]==1]+1
  # Assign n_round outside of the function

  # Decide who takes up invite
  m.Screen[ ,t+1 , "Respond_DS"] <- ((m.Rand[ ,"Respond_DS", t] < scr.Params[, "DS_uptake"]) & m.Screen[ ,t+1 , "Invite_DS"])*1

  # Decide who gets a positive result (includes true positives and false positives)
  m.Screen[ ,t+1 , "Positive_DS"] <- ((m.Rand[ ,"Positive_DS", t] < scr.Params[, "DS_diag"]) & m.Screen[ ,t+1 , "Respond_DS"]==1)*1
  
  # Decide who takes up Cystoscopy
  m.Screen[ ,t+1 , "Respond_Cyst"] <- (m.Rand[ ,"Respond_Cyst", t] < scr.Params[, "Cyst_uptake"] & m.Screen[ ,t+1 , "Positive_DS"] ==1)*1
                                       
  # Probability to be diagnosed in each state
  Cystoscopy_diag <- (m.State %*% test_accuracy[, "Sens"])
  
  # Decide who is diagnosed through cystoscopy
  m.Screen[ ,t+1 , "Diagnostic_Cyst"] <- (m.Rand[ ,"Positive_Cyst", t] < Cystoscopy_diag[,1] & m.Screen[ ,t+1 , "Respond_Cyst"]==1)*1
  
  # Decide who is diagnosed with  LG and HG tumours
  m.Screen[ ,t+1 , "LG"] <- (m.Screen[ ,t+1 , "Diagnostic_Cyst"]==1 & m.State[, "BC_LG"] ==1)*1
  m.Screen[ ,t+1 , "HG"] <- (m.Screen[ ,t+1 , "Diagnostic_Cyst"]==1 & (m.State[, "St1_HG"] ==1 | m.State[, "St2_HG"] ==1 | m.State[, "St3_HG"] ==1 | m.State[, "St4_HG"] ==1))*1
  
  #Decide who had TURBT (everyone with detected LG and HG cancer). I guess FP cystoscopy can't lead to TURBT?
  #m.Screen[ ,t+1 , "TURBT"] <- (m.Screen[ ,t+1 , "LG"]==1  | m.Screen[ ,t+1 , "HG"]==1)*1
  
  #Decide who had TURBT (everyone with detected LG and HG cancer including FP)
  m.Screen[ ,t+1 , "TURBT"] <- m.Screen[ ,t+1 , "Diagnostic_Cyst"]
  
  #Decide who died from TURBT
  m.Screen[ ,t+1 , "Die_TURBT"] <- ((m.Rand[ ,"Die_TURBT", t] < Mort.TURBT) & m.Screen[ ,t+1 , "TURBT"]==1)*1
    
  #Decide who had FP (everyone without cancer and FP dipstick and cystoscopy)
  m.Screen[ ,t+1 , "FP"] <- (m.Screen[ ,t+1 , "Diagnostic_Cyst"]==1 & m.State[, "NoBC"] ==1)*1
  
  #Decide who had FN (everyone with cancer who hasn't been diagnosed)
  m.Screen[ ,t+1 , "FN"] <- (m.Screen[ ,t+1 , "TURBT"]==0 & (m.State[, "St1_HG"] ==1 | m.State[, "St2_HG"] ==1 | m.State[, "St3_HG"] ==1 | m.State[, "St4_HG"] ==1))*1
  
  m.Screen
  
}


#' @details
#' This function updates m.Diag following screening 
#' @params
#' m:Screen: an array of screening and surveillance history
#' m.State: a matrix of current health states for all individuals
#' m.Diag: a matrix giving diagnostic status for all individuals
#' pop: a population matrix of individuals, each with a current age
#' t: current time point
#' scr.Params: updated screening parameters
#' @return an updated current diagnostic information matrix
#' 
f.screen_diag <- function(m.Screen, m.State, m.Diag, pop, t) {
  
  #Modify diagnosis columns to include new BC diagnoses from screening and surveillance
  m.Diag[, "HG_new_diag"] <- m.Diag[, "HG_new_diag"] +m.Screen[ ,t+1 , "HG"]
  m.Diag[, "HG_yr_diag"] <- m.Diag[, "HG_yr_diag"] + m.Screen[ ,t+1 , "HG"]
  m.Diag[, "LG_yr_diag"] <- m.Diag[, "LG_yr_diag"] + m.Screen[ ,t+1 , "LG"]
  m.Diag[, "yr_diag"] <- m.Diag[, "yr_diag"] + m.Screen[ ,t+1 , "HG"] + m.Screen[ ,t+1 , "LG"]
  
  m.Diag[, "HG_diag"] <- m.Diag[, "HG_diag"] + m.Screen[ ,t+1 , "HG"]
  m.Diag[, "HG_screen_diag"] <- m.Diag[, "HG_screen_diag"] + m.Screen[ ,t+1 , "HG"]
  m.Diag[, "HG_age_diag"] <- m.Diag[, "HG_age_diag"] + (pop[, "age"] * m.Screen[ ,t+1 , "HG"])
  m.Diag[, "HG_stage_diag"] <- m.Diag[, "HG_stage_diag"] + 
    ((m.State[, "St1_HG"] * 1 + m.State[, "St2_HG"] * 2 + m.State[, "St3_HG"] * 3+ m.State[, "St4_HG"] * 4)* m.Screen[ ,t+1 , "HG"])
  
  m.Diag[, "LG_new_diag"] <- m.Diag[, "LG_new_diag"] +m.Screen[ ,t+1 , "LG"]
  m.Diag[, "LG_diag"] <- m.Diag[, "LG_diag"] + m.Screen[ ,t+1 , "LG"]
  m.Diag[, "LG_screen_diag"] <- m.Diag[, "LG_screen_diag"] + m.Screen[ ,t+1 , "LG"]
  m.Diag[, "LG_age_diag"] <- m.Diag[, "LG_age_diag"] + (pop[, "age"] * m.Screen[ ,t+1 , "LG"])
  
  m.Diag
}


#' @details
#' This function updates m.M following screening 
#' @params
#' m:Screen: an array of screening and surveillance history
#' m.State: a matrix of current health states for all individuals
#' m.Diag: a matrix giving diagnostic status for all individuals
#' m.M: a matrix with health states for each individual
#' pop: a population matrix of individuals, each with a current age
#' t: current time point
#' scr.Params: updated screening parameters
#' @return an updated current m.M
#' 
#' 

#f.update.m.M <- function(m.M, m.Screen, m.Diag, t){
  
  # Replace with death for those who died from perforation during TURBT
 # m.M[,t+1][which(m.Screen[ ,t+1 , "Die_TURBT"]==1)] <- 4
  
  # chloe: there needs to be an adjustment if patients need to be assumed to return in no cancer state in the end of the survival data -10 y for HRBC
  # Replace with no cancer for those who had HG cancer but survived for 10 years
 # m.M[,t+1][which(m.Diag[ ,"HG_yr_diag"] >=10)] <- 1
  
  # Replace with no cancer for those who had LG cancer for 3 years
 # m.M[,t+1][which(m.Screen[ ,t+1 , "LG"]==1 & m.Diag[ ,"LG_yr_diag"] >=3)] <- 1
  
 # return(m.M)
#}
