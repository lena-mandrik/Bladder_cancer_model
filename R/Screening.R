### Functions for screening ###

#' @details
#' This function calculates individualised screening parameters for each person
#' @params
#' pop: a population matrix of individuals, each with a current age
#' m:Screen: an array of screening history
#' m.State: a matrix giving current health state for all individuals
#' @return a matrix of individualised screening parameters
#' 

f.calc.screen.params <- function(pop, m.Screen, m.M, t, m.State, m.State.KC, disease) {
  
    #Calculates dipstick uptake by personal characteristics
    DS_uptake <- 1/(1+exp(-((cbind(m.M[,t] !=5, pop[,"age"] <55, 
                                    pop[,"age"] >=55 & pop[,"age"] <60, pop[,"age"] >=65 & pop[,"age"] <70, pop[,"age"] >=70, 
                                    pop[, "sex"] ==0, (rowSums(m.Screen[, ,"Invite_DS"]) >=1 & rowSums(m.Screen[, ,"Respond_DS"]) ==0),
                                    (rowSums(m.Screen[, ,"Respond_DS"]) >=1), pop[, "imd"] ==2, pop[, "imd"] ==3, pop[, "imd"] ==4, 
                                    pop[, "imd"] ==5, pop[, "ethnic"] ==3) * 1) %*% coef_DT_Uptk)))
    
    colnames(DS_uptake) <- "DS_uptake"
    
    #Add diagnostic uptake: considered as 1 in the basecase analysis
    Diag_uptake <- rep(Diag.UPTK, nsample)
    
    if(disease =="bladder_kidney"){
      Diag_kidney =f.diag.params(m.State.KC, .env=e.KC, "KC")
      Diag_bladder =f.diag.params(m.State, .env=e.BC, "BC")
      scr.Params <- cbind(DS_uptake, Diag_uptake, Diag_bladder, Diag_kidney)
    } else if(disease =="bladder"){
      Diag_bladder =f.diag.params(m.State, .env=e.BC, "BC")
      scr.Params <- cbind(DS_uptake, Diag_uptake, Diag_bladder)
    } else{
      Diag_kidney =f.diag.params(m.State.KC, .env=e.KC, "KC")
      scr.Params <- cbind(DS_uptake, Diag_uptake, Diag_kidney)
    }
    
    return(scr.Params)
  }

################################################################################################################################
#' @details
#' This function calculates individualised screening parameters for each person
#' USed as an input to f.calc.screen.params 
#' @params
#' .env: environment either e.KC (kidney cancer) or e.BC (bladder cancer)
#' m.State: a matrix giving current health state for all individuals
#' @return a matrix of individualised screening parameters
#' 
f.diag.params <- function(m.State, .env, name.C){
  
  #Calculates DT sensitivity or false positives depending upon underlying health state and the modelled disease
  # for bladder cancer
  DS_diag <- (m.State %*% .env$test_accuracy[, "Sens"]) 
  
  
  # Probability to be diagnosed in each state at primary care (US + cytology)
  GP_diag <- (m.State %*% .env$diag1_accuracy[, "Sens"])
  
  # Probability to be diagnosed in each state at haematuria center (Cystoscopy)
  Cystoscopy_diag <- (m.State %*% .env$diag2_accuracy[, "Sens"])
  
  scr.Params <- cbind(DS_diag, GP_diag, Cystoscopy_diag)
  
  colnames(scr.Params) <- paste0(c("DS_diag", "GP_diag", "Cyst_diag"), "_", name.C)
  
  return(scr.Params)
  
}

################################################################################################################################


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
f.DS_screen <- function(m.Screen, m.Diag, m.M, m.Rand, pop, t, scr.Params, DS_age, DS_round, DS_freq, n_round, disease) {

  # set n_round as the current screening round
  
  # set eligibility by smoking status (considering that smoking status is variable)
  smoke_eligible <- as.matrix(rep(1, nsample), ncol=1) #set eligibility of the population by their current smoking status on that year
  if(screen_elig == "c.smoke"){
    smoke_eligible[pop[ ,"current_smoke"] !=1, ] <- 0
    }else if(screen_elig == "all.smoke"){
      smoke_eligible[pop[ ,"no_smoke"] ==1 , ] <- 0
    }
  
  if(cohort==1){
    m.Screen[ ,t+1 , "Invite_DS"] <- (m.Diag[, "LG_diag"] ==0 & m.Diag[, "Diag"] ==0 & # population not diagnosed with cancer yet
                                        m.M[,t] !=5 & # alive population
                                        smoke_eligible[,1]==1 & # the pop status by their current smoking eligibility
                                        (pop[, "age"] == DS_age | # for one-time screening, should be equal to the pop age; or the screening start for the repetitive screening
                                        (pop[, "age"] > DS_age & n_round < DS_round & (DS_freq ==1 | # for annual screening
                                        (m.Screen[ ,t, "Invite_DS"]==0 & DS_freq ==2)))))*1 # for biennial screening
  } else if(cohort ==0){
    
    # Not updated with the ch=ohort ==0
    m.Screen[ ,t+1 , "Invite_DS"] <- (m.Diag[, "LG_diag"] ==0 & m.Diag[, "HG_diag"] ==0 & # population not diagnosed with cancer yet
                                        m.State[, "DeathBC"] ==0 & m.State[, "DeathOC"] ==0 & # alive population
                                        smoke_eligible[,1]==1 & #the pop status by their current smoking eligibility
                                        (t==1 | # for one-time screening, should be equal to the pop age; or the screening start for the repetitive screening
                                           (n_round < DS_round & (DS_freq ==1 | # for annual screening
                                            (m.Screen[ ,t, "Invite_DS"]==0 & DS_freq ==2)))))*1 # for biennial screening
  }

  
  # Update the number of screen rounds each person received 
  n_round[m.Screen[ ,t+1 , "Invite_DS"]==1] <- n_round[m.Screen[ ,t+1 , "Invite_DS"]==1]+1
  # Assign n_round outside of the function

  # Decide who takes up invite
  m.Screen[ ,t+1 , "Respond_DS"] <- ((m.Rand[ ,"Respond_DS", t] < scr.Params[, "DS_uptake"]) & m.Screen[ ,t+1 , "Invite_DS"])*1

  # Decide who gets a positive result (includes true positives and false positives)
  m.Screen[ ,t+1 , "Positive_DS"] <- ((m.Rand[ ,"Positive_DS", t] < scr.Params[, "DS_diag_BC"] | m.Rand[ ,"Positive_DS", t] < scr.Params[, "DS_diag_KC"]) & m.Screen[ ,t+1 , "Respond_DS"]==1)*1
  
  # Decide who takes up invite to the follow-up of DS positive test
  # Set up for now as one uptake for diagnostic similar to cystoscopy 
  m.Screen[ ,t+1 , "Respond_diag"] <- ((m.Rand[ ,"Respond_diag", t] < scr.Params[, "Diag_uptake"]) & m.Screen[ ,t+1 , "Positive_DS"]==1)*1
  
  # Decide who gets a positive result at GP follow up visit (US  + cytology)
  m.Screen[ ,t+1 , "Positive_diag"] <- ((m.Rand[ ,"Positive_diag", t] < scr.Params[, "GP_diag_BC"] | m.Rand[ ,"Positive_diag", t] < scr.Params[, "GP_diag_KC"]) & m.Screen[ ,t+1 , "Respond_diag"]==1)*1
  
  # Decide who takes up Cystoscopy
  m.Screen[ ,t+1 , "Respond_Cyst"] <- (m.Rand[ ,"Respond_Cyst", t] < scr.Params[, "Diag_uptake"] & m.Screen[ ,t+1 , "Positive_diag"] ==1)*1
                                       
  # Decide who is diagnosed through cystoscopy
  m.Screen[ ,t+1 , "Diagnostic_Cyst"] <- ((m.Rand[ ,"Positive_Cyst", t] < scr.Params[, "Cyst_diag_BC"] | m.Rand[ ,"Positive_Cyst", t] < scr.Params[, "Cyst_diag_KC"]) & m.Screen[ ,t+1 , "Respond_Cyst"]==1)*1
  
  # Decide who is diagnosed with  LG, HG, KC tumours
  m.Screen[ ,t+1 , "LG"] <- (m.Screen[ ,t+1 , "Diagnostic_Cyst"]==1 & m.M[,t+1]==2)*1
  m.Screen[ ,t+1 , "HG"] <- (m.Screen[ ,t+1 , "Diagnostic_Cyst"]==1 & m.M[,t+1]==3)*1
  m.Screen[ ,t+1 , "KC"] <- (m.Screen[ ,t+1 , "Diagnostic_Cyst"]==1 & m.M[,t+1]==4)*1
  
  #Decide who had TURBT (everyone with detected LG and HG cancer including FP)
  m.Screen[ ,t+1 , "Surgery"] <- m.Screen[ ,t+1 , "Diagnostic_Cyst"]
  
  # Assign a probability to die during surgery based on which disease is modelled
  if(disease=="bladder" | disease=="bladder_kidney"){Mort.TURBT=e.BC$Mort.TURBT} else{Mort.TURBT=e.KC$Mort.TURBT}
  
  #Decide who died from TURBT
  m.Screen[ ,t+1 , "Die_Surgery"] <- ((m.Rand[ ,"Die_Surgery", t] < Mort.TURBT) & m.Screen[ ,t+1 , "Surgery"]==1)*1
    
  #Decide who had FP (everyone without cancer and FP dipstick and cystoscopy)
  m.Screen[ ,t+1 , "FP_BC"] <- (m.Screen[ ,t+1 , "Diagnostic_Cyst"]==1 & m.Rand[ ,"Positive_Cyst", t] < scr.Params[, "Cyst_diag_BC"] & m.M[,t+1]==1)*1
  m.Screen[ ,t+1 , "FP_KC"] <- (m.Screen[ ,t+1 , "Diagnostic_Cyst"]==1 & m.Rand[ ,"Positive_Cyst", t] < scr.Params[, "Cyst_diag_KC"] & m.M[,t+1]==1)*1
  
  #Decide who had FN (everyone with cancer who hasn't been diagnosed)
  m.Screen[ ,t+1 , "FN"] <- (m.Screen[ ,t+1 , "Diagnostic_Cyst"]==0 & (m.M[,t+1]==2 | m.M[,t+1]==3 | m.M[,t+1]==4))*1
  
  
  
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
f.screen_diag <- function(m.Screen, m.State, m.State.KC, m.Diag, pop, t) {
  
  #Modify diagnosis columns to include new cancer diagnoses from screening and surveillance
  yr_diag <- m.Diag[, "yr_diag"] + m.Screen[ ,t+1, "HG"]+ m.Screen[ ,t+1 , "KC"]
  
  # Consider FP as LG, as it is more likely to be LG than HG
  m.Diag[, "LG_yr_diag"] <- m.Diag[, "LG_yr_diag"] + m.Screen[ ,t+1, "LG"] #+ m.Screen[ ,t+1 , "FP"]
  
  # Record FP for both kidney and bladder cancers
  m.Diag[, "FP_BC"] <- m.Diag[, "FP_BC"] + m.Screen[ ,t+1, "FP_BC"]
  m.Diag[, "FP_KC"] <- m.Diag[, "FP_KC"] + m.Screen[ ,t+1, "FP_KC"]
  
  # Record diagnostic parameters fro those diagnosed through screening
  m.Diag[, "Diag"] <- m.Diag[, "Diag"] + m.Screen[ ,t+1 , "HG"] + m.Screen[ ,t+1 , "KC"]
  m.Diag[, "Screen_diag"] <- m.Diag[, "Screen_diag"] + m.Screen[ ,t+1 , "HG"] + m.Screen[ ,t+1 , "KC"]
  m.Diag[, "Age_diag"] <- m.Diag[, "Age_diag"] + (pop[, "age"] * m.Screen[ ,t+1 , "HG"])+ (pop[, "age"] * m.Screen[ ,t+1 , "KC"])
  
  m.Diag[, "Stage_diag"] <- m.Diag[, "Stage_diag"] + 
    ((m.State[, "St1_HG"] * 1 + m.State[, "St2_HG"] * 2 + m.State[, "St3_HG"] * 3+ m.State[, "St4_HG"] * 4)* m.Screen[ ,t+1 , "HG"])+
  ((m.State.KC[, "St1_HG"] * 1 + m.State.KC[, "St2_HG"] * 2 + m.State.KC[, "St3_HG"] * 3+ m.State.KC[, "St4_HG"] * 4)* m.Screen[ ,t+1 , "KC"])
  
  m.Diag[, "LG_diag"] <- m.Diag[, "LG_diag"] + m.Screen[ ,t+1 , "LG"] #+m.Screen[ ,t+1 , "FP_BC"]
  m.Diag[, "LG_screen_diag"] <- m.Diag[, "LG_screen_diag"] + m.Screen[ ,t+1 , "LG"]#+m.Screen[ ,t+1 , "FP"]
  m.Diag[, "LG_age_diag"] <- m.Diag[, "LG_age_diag"] + (pop[, "age"] * m.Screen[ ,t+1 , "LG"]) #+ (pop[, "age"] *m.Screen[ ,t+1 , "FP"])
  
  m.Diag[, "yr_diag"] <- m.Diag[, "LG_yr_diag"]
  m.Diag[, "yr_diag"] <- replace(m.Diag[, "yr_diag"], yr_diag >0, yr_diag[yr_diag>0])
  
  m.Diag
}



#######################################################
#re-set the stage at m.M_8s if the patient stage changed at this cycle and progression was to happen in the 2d half of the year
f.screen.shift <- function(m.M_8s, m.C.T.to.Stage, m.Screen, m.Diag, t, c.name){
  
  # concervative way of downstaging one step down
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.M_8s[, t+1]==5 & ceiling(m.C.T.to.Stage[ ,"T.onsetToStage2"])==m.Diag[,"yr_onset"] 
                           & m.C.T.to.Stage[ ,"T.onsetToStage2"]%%1 > 0.5 & m.Screen[ ,t+1 , c.name]==1, 3) 
  #get downstage to one stage (previously sampled stage) if time to progression is less than 6 m
  
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.M_8s[, t+1]==6 & ceiling(m.C.T.to.Stage[ ,"T.onsetToStage3"])==m.Diag[,"yr_onset"] 
                           & m.C.T.to.Stage[ ,"T.onsetToStage3"]%%1 > 0.5 & m.Screen[ ,t+1 , c.name]==1, 4) 
  
  #get downstage (previously sampled stage) if time to progression is less than 6 m
  m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.M_8s[, t+1]==7 & ceiling(m.C.T.to.Stage[ ,"T.onsetToStage4"])==m.Diag[,"yr_onset"] 
                           & m.C.T.to.Stage[ ,"T.onsetToStage4"]%%1 > 0.5 & m.Screen[ ,t+1 , c.name]==1, 6) 
  # get downstage (previously sampled stage) if time to progression is less than 6 m

 
  return(m.M_8s)
  }