### Function to run the simulation ###

Simulate_NHD <- function(nsample, n.t, pop, m.Rand) { 
  
  # Arguments:     
  # nsample:     number of individuals simulated; Global parameters, Master script
  # n.t:     total number of cycles to run the model; Global parameters, Master script   
  # pop:     matrix of individual level population attributes
  # Makes use of:   
  # Probs:  function for the estimation of transition probabilities (script "Functions")
  # Costs:   function for the estimation of cost state values (script "Functions")
  # Effs:    function for the estimation of state specific health outcomes (QALYs) (script "Functions")   
  # create matrices capturing the state name/costs/health outcomes for all individuals at each time point 
  
  model_matrices = initialize_matrices(nsample, n.t, n.s_long, states_long, run_mode)
  m.M = model_matrices$m.M
  m.M_8s = model_matrices$m.M_8s
  m.M_8s_KC =model_matrices$m.M_8s_KC
  m.C =model_matrices$m.C
  m.E =model_matrices$m.E
  m.Diag=model_matrices$m.Diag
  m.State =model_matrices$m.State
  m.State.KC =model_matrices$m.State.KC
  m.Screen =model_matrices$m.Screen
  m.Out=model_matrices$m.Out
  
  # Set the number of screening rounds each person receives; counts as zero before the cycles
  n_round <<- rep(0, nsample) 
  
  # loop to run the model over time
  for(t in 1:n.t) {
    
    # Natural History
    TP <- f.calc.indiv.TPs(pop, m.Diag, disease) #calculate new individualised transition probabilities for onset of BC and OC mortality
    
    #if(t==1){
     # TP[c(1:5),3]=1
      #TP[c(6:15),5]=1
     # TP[c(1:5),c(1,2,4,5)]=0
     # TP[c(6:15),c(1:4)]=0
     # TP[c(16:20),2]=1
     # TP[c(16:20),c(1,3:5)]=0 
  #  }
   
    
    m.p <- f.Probs(m.M[, t], TP) #calculate transition probabilities for 4 states at cycle t (excludes TP for those with invasive BC)
    
    m.M[, t+1] <- f.samplev(m.p, m.Rand, t) #Sample the next health state and store it in the Matrix m.M numerically
    # 1- no cancer, 2- LGBC cancer, 3 -HGBC cancer, 4 -KC, 5- mortality
    
    # update the m.Diag matrix with the LR and HR BC or KC onset 
    m.Diag <- f.onset.Diag(m.Diag, m.M, pop, t) #includes m.Diag[, "HG_age_onset"]
    
    # Update m.M_8s matrix with the states for those who have HGBC or KC
    if(disease=="bladder" | disease=="bladder_kidney"){
      m.M_8s <- f.C.stages(m.M_8s, m.BC.T.to.Stage, m.M, m.Diag, t, "bladder")
        
    } 
    if(disease=="kidney" | disease=="bladder_kidney"){
      m.M_8s_KC <- f.C.stages(m.M_8s_KC, m.KC.T.to.Stage, m.M, m.Diag, t, "kidney")
    }

    # Update the matrix with the health state numbers with either 0 or 1 depending if it is equal to the sampled state
    m.State[] <- m.State.KC[] <- 0
    
    for(n in 1:n.s_long) {
      if(disease=="bladder" | disease=="bladder_kidney"){m.State[, n] <- replace(m.State[, n], m.M_8s[, t+1] ==n, 1)}
      if(disease=="kidney" | disease=="bladder_kidney"){m.State.KC[, n] <- replace(m.State.KC[, n], m.M_8s_KC[, t+1] ==n, 1)}
    }
    
    # Symptomatic detection
    # If run with screening (annual cycles) consider that 50% of patients become symptomatic the 1st half of the year if they are symptomatic that year
    if(DS_screen ==1){elig_time <- as.matrix(1*(m.Rand[,"Screen_time", t] <= 0.5), ncol=1)} else {elig_time <- as.matrix(rep(1,nsample), ncol=1)}  # Run symptomatic mode for those with m.Rand < 0.5 (half of the people) to ensure equal impact of screen and sympt pathways
    
    # move year of diagnosis on by one
    m.Diag[m.M[,t] != 5, "yr_diag"][m.Diag[m.M[,t] != 5, "yr_diag"] >= 1] <- m.Diag[m.M[,t] != 5, "yr_diag"][m.Diag[m.M[,t] != 5, "yr_diag"] >= 1] + 1
    m.Diag[m.M[,t] != 5, "LG_yr_diag"][m.Diag[m.M[,t] != 5, "LG_yr_diag"] >= 1] <- m.Diag[m.M[,t] != 5, "LG_yr_diag"][m.Diag[m.M[,t] != 5, "LG_yr_diag"] >= 1] + 1
    
    m.Diag <- f.symptom(m.Diag, m.State, m.State.KC, m.Rand, pop, t, m.M, elig_time, nsample) 
  
    # Screening detection
    if(DS_screen ==1){
      
      scr.Params <- f.calc.screen.params(pop, m.Screen, m.M, t, m.State, m.State.KC, disease) 

      m.Screen <- f.DS_screen(m.Screen, m.Diag, m.M, m.Rand, pop, t, scr.Params, DS_age, DS_round, DS_freq, n_round, disease)
      
      #re-set the stage at m.M_8s if the patient stage changed at this cycle and progression was to happen in the 2d half of the year
      if(disease=="bladder" | disease=="bladder_kidney"){m.M_8s <- f.screen.shift(m.M_8s, m.BC.T.to.Stage, m.Screen, m.Diag ,t, "HG")}
      if(disease=="kidney" | disease=="bladder_kidney"){m.M_8s_KC <- f.screen.shift(m.M_8s_KC, m.KC.T.to.Stage, m.Screen, m.Diag ,t, "KC")}
      
      # update again state matrix
      m.State[] <- m.State.KC[] <- 0
      for(n in 1:n.s_long) {
        if(disease=="bladder" | disease=="bladder_kidney"){m.State[, n] <- replace(m.State[, n], m.M_8s[, t+1] ==n, 1)}
        if(disease=="kidney" | disease=="bladder_kidney"){m.State.KC[, n] <- replace(m.State.KC[, n], m.M_8s_KC[, t+1] ==n, 1)}
      }
     
      m.Diag <- f.screen_diag(m.Screen, m.State, m.State.KC, m.Diag, pop, t)
      
      # Replace with death for those who died from perforation during TURBT
      m.M[,t+1][which(m.Screen[ ,t+1 , "Die_Surgery"]==1)] <- 5
      
      # NOTE: currently only a progress to HGBC after the surveillance for 3 years is considered for LGBC. HGBC are only following the survival curve.
      # Not sure whether I need to return to no cancer everyone in the end of the 10 year time period

      elig_time <- as.matrix(1*(m.Rand[,"Screen_time", t] > 0.5), ncol=1)
      m.Diag <- f.symptom(m.Diag, m.State, m.State.KC, m.Rand, pop, t, m.M, elig_time, nsample) 
    }
    
    # Define deaths from cancer
    died_undiagnosed.BC <- Death_all.BC <- died_undiagnosed.KC <- Death_all.KC <- rep(0, nsample)
      
    # define BC deaths
    if(disease=="bladder" | disease=="bladder_kidney"){
    death.BC <- f.return.C.death(pop, m.Diag, m.State, e.BC, m.Rand, t)
    died_undiagnosed.BC <- death.BC$died_undiagnosed
    Death_all.BC <- death.BC$C_death_all
    }
    
    if(disease=="kidney" | disease=="bladder_kidney"){
    # define KC deaths
    death.KC <- f.return.C.death(pop, m.Diag, m.State.KC, e.KC, m.Rand, t)
    died_undiagnosed.KC <- death.KC$died_undiagnosed
    Death_all.KC <- death.KC$C_death_all
    }
    
    # update with the age of BC death
    # Record the age of death for those with cancer
    m.Diag[, "age_C_death"] <- m.Diag[, "age_C_death"] + (pop[, "age"] * Death_all.BC)  + (pop[, "age"] * Death_all.KC) 
    
    # If died undiagnosed, clear all characteristics
    m.Diag[(died_undiagnosed.BC[] ==1 | died_undiagnosed.KC[] ==1), "Diag"] <- m.Diag[(died_undiagnosed.BC[] ==1 | died_undiagnosed.KC[] ==1 ), "Sympt_diag"] <- 
      m.Diag[(died_undiagnosed.BC[] ==1 | died_undiagnosed.KC[] ==1), "yr_diag"] <-  m.Diag[(died_undiagnosed.BC[] ==1 | died_undiagnosed.KC[] ==1), "Age_diag"] <-0 

    # Update the mortality for BC (replace with the state 8 in the m.M_8s matrix and state 4 in m.M matrix those who died with BC)
    # In the rare even of the BC perforation during screening, if this patient meant to die from BC during the same cycle it would be counted as BC death
    m.M_8s[, t+1] <- ifelse(died_undiagnosed.BC > 0, 5, ifelse(Death_all.BC > 0, 8, m.M_8s[, t+1])) # If died from BC replace to state 8, if died undiagnosed - replace with 4 (OCM)
    m.M_8s_KC[, t+1] <- ifelse(died_undiagnosed.KC > 0, 5, ifelse(Death_all.KC > 0, 8, m.M_8s_KC[, t+1])) # If died from KC replace to state 8, if died undiagnosed - replace with 4 (OCM)
    
    # Update m.M matrix with deaths
    m.M[, t+1] <- replace(m.M[, t+1], m.M_8s[, t+1]==8 | m.M_8s_KC[, t+1]==8 | died_undiagnosed.BC >0 | died_undiagnosed.KC >0, 5)
    
    # Update the output matrices
    if(run_mode != "Calibration" ){
    # Update the QALYs and costs (not half cycle corrected)
    m.E[, t+1] <- f.calc.utility(m.State, m.State.KC, m.Diag, pop, t, m.M, disease) # Assess effects per individual during the cycle t+1 
    m.C[, t+1] <- f.calc.cost(m.State, m.State.KC, m.Diag, disease) # Assess CRC treatment costs per individual during the cycle t+1 (note not half cycle corrected) 

    # Gather outcomes for cycle t (total and subgroups)
    m.Out <- f.aggregate.outcomes(m.M_8s, m.M_8s_KC, m.Out, m.M, m.E, m.C, m.Diag, m.Screen, pop, t, disease)
    }
    
    # Update the age for only those individuals who are alive
    IND <- m.M[, t+1] != 5
    pop[IND, "age"] <- pop[IND,"age"] +1 # update the age
    
    ncol_pop <- ncol(pop)
    rand.quit <- m.Rand[ ,"Smoke_quit", t]
    
    if(!is.null(dim(pop[IND, ]))) {
    pop[IND, ] <- f.smoke.change(pop[IND, ], rand.quit[IND])  # update smoking status
    
    if(disease=="bladder_kidney" | disease=="bladder"){
      pop[IND, ] <- f.risk.calc(pop[IND, ], e.BC, "bladder") #update the risk of BC and p of onset of BC
    } 
     if(disease=="bladder_kidney" | disease=="kidney"){
       pop[IND, ] <- f.risk.calc(pop[IND, ], e.KC, "kidney")
     }
     
    }
    
  }#this is a loop for the time point
  

  # Extract the matrices for technical validity
  if(run_mode == "Testing") { # Create a matrix of transitions across states
    
    testing.outputs <- f.validity.matrix(m.M_8s, m.M_8s_KC, m.M, disease)
    
  } 
  
  # Extract the matrices by sex for calibration
  
  if (run_mode == "Calibration" | run_mode == "Testing") { #  add additional two matrices to report TR matrix by gender

    
    if(disease =="bladder_kidney"){
      m.sex.validity_BC <- f.validity.sex.matrix(m.M_8s); m.sex.validity_KC <- f.validity.sex.matrix(m.M_8s_KC)
      m.sex.validity <- list(TR_m.BC=m.sex.validity_BC$TR_m,
                             TR_f.BC=m.sex.validity_BC$TR_f,
                             TR_m.KC=m.sex.validity_KC$TR_m,
                             TR_f.KC=m.sex.validity_KC$TR_f)
    } else if(disease =="bladder"){
      m.sex.validity_BC <- f.validity.sex.matrix(m.M_8s)
      m.sex.validity <- list(TR_m.BC=m.sex.validity_BC$TR_m,
                             TR_f.BC=m.sex.validity_BC$TR_f)
    } else{
      m.sex.validity_KC <- f.validity.sex.matrix(m.M_8s_KC)
      m.sex.validity <- list(TR_m.KC=m.sex.validity_KC$TR_m,
                             TR_f.KC=m.sex.validity_KC$TR_f)
    }

 }
  
  
  results <- switch(run_mode, 
                    Testing = list(testing.outputs=testing.outputs, m.sex.validity=m.sex.validity, m.Diag = m.Diag, m.M_8s=m.M_8s, m.M= m.M),
                    Calibration = list(m.sex.validity=m.sex.validity, m.Diag = m.Diag),
                    Deterministic = m.Out,
                    PSA= m.Out
                    )
  
  
  return(results)  # return the results 
  }  