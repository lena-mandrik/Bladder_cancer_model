### Function to run the simulation ###

Simulate_NHD <- function(nsample, n.t, pop, m.BC.T.to.Stage) { 
  
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
  m.Screen =model_matrices$m.Screen
  m.Out=model_matrices$m.Out
  
  # Set the number of screening rounds each person receives; counts as zero before the cycles
  n_round <<- rep(0, nsample) 
  
  # loop to run the model over time
  for(t in 1:n.t) {
    
    # Natural History
    TP <- f.calc.indiv.TPs(pop, m.Diag, disease) #calculate new individualised transition probabilities for onset of BC and OC mortality
    
    TP[c(1:5),3]=1
    TP[c(6:10),5]=1
    TP[c(1:5),c(1,2,4,5)]=0
    TP[c(6:10),c(1:4)]=0
    TP[c(11:15),2]=1
    TP[c(11:15),c(1,3:5)]=0
    
    m.p <- f.Probs(m.M[, t], TP) #calculate transition probabilities for 4 states at cycle t (excludes TP for those with invasive BC)
    
    m.M[, t+1] <- f.samplev(m.p, m.Rand, t) #Sample the next health state and store it in the Matrix m.M numerically
    # 1- no cancer, 2- LGBC cancer, 3 -HGBC cancer, 4 -KC, 5- mortality
    
    # update the m.Diag matrix with the LR and HR BC or KC onset 
    m.Diag <- f.onset.Diag(m.Diag, m.M, pop, t) #includes m.Diag[, "HG_age_onset"]
    
    # Update m.M_8s matrix with the states for those who have HGBC or KC
    if(disease=="bladder" | disease=="bladder_kidney"){
      m.M_8s <- f.HG.stages(m.M_8s, m.BC.T.to.Stage, m.M, m.Diag, t)
    } 
    if(disease=="kidney" | disease=="bladder_kidney"){
      m.M_8s_KC <- f.KC.stages(m.M_8s, m.KC.T.to.Stage, m.M, m.Diag, t)
    }

    # Update the matrix with the health state numbers with either 0 or 1 depending if it is equal to the sampled state
    m.State[] <- 0
    for(n in 1:n.s_long) {
      m.State[, n] <- replace(m.State[, n], m.M_8s[, t+1] ==n, 1)
    }
    
    #Symptomatic detection
    if(DS_screen ==1){elig_time <- as.matrix(1*(m.Rand[,"Screen_time", t] <= 0.5), ncol=1)} else {elig_time <- as.matrix(rep(1,nsample), ncol=1)}  # Run symptomatic mode for those with m.Rand < 0.5 (half of the people) to ensure equal impact of screen and sympt pathways
    
    # move year of diagnosis on by one
    m.Diag[m.M[,t] != 4, "HG_yr_diag"][m.Diag[m.M[,t] != 4, "HG_yr_diag"] >= 1] <- m.Diag[m.M[,t] != 4, "HG_yr_diag"][m.Diag[m.M[,t] != 4, "HG_yr_diag"] >= 1] + 1
    m.Diag[m.M[,t] != 4, "LG_yr_diag"][m.Diag[m.M[,t] != 4, "LG_yr_diag"] >= 1] <- m.Diag[m.M[,t] != 4, "LG_yr_diag"][m.Diag[m.M[,t] != 4, "LG_yr_diag"] >= 1] + 1
    
    m.Diag <- f.symptom(m.Diag, m.State, m.Rand, pop, t, m.M, elig_time, nsample) 
    
  
    # Screening detection
    if(DS_screen ==1){
      
      scr.Params <- f.calc.screen.Params(pop, m.Screen, m.State, test_accuracy, diag1_accuracy, diag2_accuracy)
      m.Screen <- f.DS_screen(m.Screen, m.Diag, m.State, m.Rand, pop, t, scr.Params, DS_age, DS_round, DS_freq, n_round)
      
      #re-set the stage at m.M_8s if the patient stage changed at this cycle and progression was to happen in the 2d half of the year
      m.M_8s <- f.screen.shift(m.M_8s, m.BC.T.to.Stage, m.Screen, m.Diag ,t)
      
      #update again state matrix
      m.State[] <- 0
      for(n in 1:n.s_long) {
        m.State[, n] <- replace(m.State[, n], m.M_8s[, t+1] ==n, 1)
      }
      
      m.Diag <- f.screen_diag(m.Screen, m.State, m.Diag, pop, t)
      
      # Replace with death for those who died from perforation during TURBT
      m.M[,t+1][which(m.Screen[ ,t+1 , "Die_TURBT"]==1)] <- 4
      # NOTE: currently only a progress to HGBC after the surveillance for 3 years is considered for LGBC. HGBC are only following the survival curve.
      # Not sure whether I need to return to no cancer everyone in the end of the 10 year time period

      elig_time <- as.matrix(1*(m.Rand[,"Screen_time", t] > 0.5), ncol=1)
      m.Diag <- f.symptom(m.Diag, m.State, m.Rand, pop, t, m.M, elig_time, nsample)
    }
    
    # define BC deaths
    TP_BC.mort <- f.calc.BCmort.TP(pop, m.Diag)
    new_BC1_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.1.mort"]) & m.State[ ,"St1_HG"]>0)
    new_BC2_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.2.mort"]) & m.State[ ,"St2_HG"]>0)
    new_BC3_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.3.mort"]) & m.State[ ,"St3_HG"]>0)
    new_BC4_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.4.mort"]) & m.State[ ,"St4_HG"]>0)
    
    # TP.ungiag.dead - should be betwen -0.001 and -0.05
    # Update the mortality with a possibility for some BC at stage 4 to be undiagnosed death 
    died_undiagnosed <- 1*(m.Rand[ ,"Death_BC_undiag", t] < (1-exp(P.ungiag.dead * ((pop[, "age"] - 80) * (1*(pop[, "age"] > 80)))))
                           & m.Diag[ ,"yr_diag"]==1
                           & new_BC4_death[]==1)
    
    new_BC4_death[died_undiagnosed[] ==1] <-0    
    
    BC_death_all <- rowSums(cbind(new_BC1_death, new_BC2_death, new_BC3_death, new_BC4_death))
    
    # update with the age of BC death
    m.Diag[, "age_BC_death"] <- m.Diag[, "age_BC_death"] + (pop[, "age"] * BC_death_all) # Record the age of death for those with cancer
    m.Diag[died_undiagnosed[] ==1, "HG_diag"] <- m.Diag[died_undiagnosed[] ==1, "HG_sympt_diag"] <- 
      m.Diag[died_undiagnosed[] ==1, "HG_yr_diag"] <-  m.Diag[died_undiagnosed[] ==1, "HG_age_diag"] <-0 # Record the age of death for those with cancer

    # Update the mortality for BC (replace with the state 8 in the m.M_8s matrix and state 4 in m.M matrix those who died with BC)
    # In the rare even of the BC perforation during screening, if this patient meant to die from BC during the same cycle it would be counted as BC death
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], BC_death_all>0, 8) # If died from BC replace to state 8
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], died_undiagnosed>0, 4) # If died undiagnosed replace with 4 (OCM)
    m.M[, t+1] <- replace(m.M[, t+1], m.M_8s[, t+1]==8 | died_undiagnosed>0, 4)
    
   #print(which(died_undiagnosed >0))
    
    if(run_mode != "Calibration" ){
    # Update the QALYs and costs (not half cycle corrected)
    m.E[, t+1] <- f.calc.utility(m.State, m.Diag, pop, t) # Assess effects per individual during the cycle t+1 
    m.C[, t+1] <- f.calc.cost(m.State, m.Diag, m.Cost.treat) # Assess CRC treatment costs per individual during the cycle t+1 (note not half cycle corrected) 
    
    # Gather outcomes for cycle t (total and subgroups)
    m.Out <- f.aggregate.outcomes(m.M_8s, m.Out, m.M, m.E, m.C, m.Diag, m.Screen, m.State, pop, t)
    }
    
    # Update the age for only those individuals who are alive
    IND <- m.M[, t+1] != 4
    pop[IND, "age"] <- pop[IND,"age"] +1 # update the age
    
    ncol_pop <- ncol(pop)
    rand.quit <- m.Rand[ ,"Smoke_quit", t]
    
    if(!is.null(dim(pop[IND, ]))) {
    pop[IND, ] <- f.smoke.change(pop[IND, ], rand.quit[IND])  # update smoking status
    pop[IND, ] <- f.risk.calc(pop[IND, ]) #update the risk of BC and p of onset of BC
    
    }
    
  }#this is a loop for the time point
  

  # Extract the matrices for technical validity
  if(run_mode == "Testing") { # Create a matrix of transitions across states
    TS_8 <- paste(m.M_8s, cbind(m.M_8s[, -1], NA), sep = "->") # transitions from one state to the other     
    TS_8 <- matrix(TS_8, nrow = nsample)
    TS_4 <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other     
    TS_4 <- matrix(TS_4, nrow = nsample)
    rownames(TS_8) <- rownames(TS_4) <-paste("Ind",   1:nsample, sep = " ")   # name the rows      
    colnames(TS_8) <- colnames(TS_4) <-paste("Cycle", 0:n.t, sep = " ")   # name the columns    
  
    TR_8 <- t(apply(m.M_8s, 2, function(x) table(factor(x, levels = v.n_long, ordered = TRUE))))     
    TR_8 <- TR_8 / nsample    # create a distribution trace
    TR_4 <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))     
    TR_4 <- TR_4 / nsample    # create a distribution trace
    
    rownames(TR_8) <- rownames(TR_4) <-paste("Cycle", 0:n.t, sep = " ")  # name the rows
    colnames(TR_8) <- v.n_long     # name the columns    
    colnames(TR_4) <-v.n    # name the columns    
    
  } else {   
    TS_8 <- TS_4 <- TR_8 <- TR_4 <- NULL   
  } 
  
  # Extract the matrices by sex for calibration
  
  if (run_mode == "Calibration" | run_mode == "Testing") { #  add additional two matrices to report TR matrix by gender
    names_Diag <- c(colnames(m.Diag), "sex", "age")
    m.Diag <- cbind(m.Diag, pop[ ,"sex"], pop[ ,"age_0"])
    colnames(m.Diag) <-names_Diag
    
    a.Temp <- array(data = NA, dim = c(nsample, n.t+1, n.s_long))
    for (i in 1:n.s_long) {
      a.Temp[ , ,i] <- (m.M_8s == i)
    }
    
    if(cohort ==1){
      a.Temp <- a.Temp*pop[,"weighting"]
    }
    
    a.Male <- a.Temp[pop[, "sex"] == 1, ,]
    if(cohort ==1){
      TR_m <- as.matrix(colSums(a.Male))/ sum(pop[, "weighting"]*pop[, "sex"])
    } else {TR_m <- as.matrix(colSums(a.Male))/ sum(pop[, "sex"])}
    rownames(TR_m) <- paste("Cycle", 0:n.t, sep = " ")  # name the rows
    colnames(TR_m) <- v.n_long     # name the columns    
    
    a.Female <- a.Temp[pop[ ,"sex"] == 0, ,]
    if(cohort ==1){
      TR_f <- as.matrix(colSums(a.Female))/ sum(pop[, "weighting"]*(1-pop[, "sex"]))
    } else {TR_f <- as.matrix(colSums(a.Female))/ sum((1-pop[, "sex"]))}
    rownames(TR_f) <- paste("Cycle", 0:n.t, sep = " ")  # name the rows
    colnames(TR_f) <- v.n_long     # name the columns 
    
  }
  
  if(run_mode == "Testing") {
    cycle <- pop[ ,"age_0"] -30
    m.M_8s <- cbind(pop[, "PID"],pop[, "age_0"], cycle, m.M_8s)
    
  }
  
  results <- switch(run_mode, 
                    Testing = list(TR_8 = TR_8, TR_4 = TR_4,TR_f = TR_f, TR_m=TR_m, m.Diag = m.Diag, m.M_8s=m.M_8s, m.M= m.M),
                    Calibration = list(TR_f = TR_f, TR_m=TR_m, m.Diag = m.Diag),
                    Deterministic = m.Out,
                    PSA= m.Out
                    )
  
  
  return(results)  # return the results 
  }  