### Function to run the simulation ###

Simulate_NHD <- function(n.i, n.t, pop) { 
  
  # Arguments:     
  # n.i:     number of individuals; Global parameters, Master script
  # n.t:     total number of cycles to run the model; Global parameters, Master script   
  # pop:     matrix of individual level population attributes
  # Makes use of:   
  # Probs:  function for the estimation of transition probabilities (script "Functions")
  # Costs:   function for the estimation of cost state values (script "Functions")
  # Effs:    function for the estimation of state specific health outcomes (QALYs) (script "Functions")   
  
  # create matrices capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.M_8s <-m.C <- m.E <-  matrix(nrow = n.i, ncol = n.t + 1)
  
  m.M[, 1] <- m.M_8s[, 1] <- rep(1, n.i)   # indicate the initial health state
  m.C[, 1] <- 0 # estimate costs per individual for the initial health state 
  m.E[, 1] <- pop[, "EQ5D"] - (Utility.age * (pop[, "age"] - pop[, "age_0"])) # estimate QALY per individual for the initial health state and starting age
  m.E[, 1] <- replace(m.E[, 1], m.E[, 1] >1, 1) # if EQ-5D over 1 reset to 1
  m.E[, 1] <- replace(m.E[, 1], m.E[, 1] <=(-0.594), -0.594) # if EQ-5D under -0.594 reset to -0.594
  
  # create another matrix capturing the current health state for each individual
  m.State <- matrix(0, nrow = n.i, ncol = n.s_long)
  colnames(m.State) <- states_long
  for(n in 1:n.s) {
    m.State[, n] <- replace(m.State[, n], m.M[, 1] ==n, 1)
  }
  
  #Create another matrix for current diagnostic information
  m.Diag <- matrix(0, nrow = n.i, ncol = 20)
  
  # When BC said, it means HG

 colnames(m.Diag) <- c("HG_state", "LG_state", "HG_diag", "LG_diag", "HG_age_diag", "LG_age_diag", "HG_sympt_diag", "LG_sympt_diag", "HG_screen_diag", "LG_screen_diag", "HG_new_diag", "LG_new_diag",
   "HG_stage_diag", "yr_diag", "HG_yr_diag", "LG_yr_diag", "HG_yr_onset", "HG_age_onset", "LG_age_onset", "age_BC_death")
  
  
  # Create an array to gather screening and surveillance information for each cycle
  m.Screen <- array(data = 0, dim = c(n.i, n.t + 1, length(screen_names)))
  dimnames(m.Screen)[[3]] <- screen_names # note that this is defined in set_params
  
  #Create a matrix for gathering aggregate outcomes
  m.Out <- matrix(0, nrow = length(out_names), ncol = n.t + 1)
  rownames(m.Out) <- out_names
  colnames(m.Out) <- c(0:n.t)
  
  #Create additional matrices for subgroup results
  #m.Out_M <- m.Out_F <- m.Out_smoke <- m.Out_past.smoke <- m.Out
  
  # loop to run the model over time
  for(t in 1:n.t) {
    
    #Natural History
    TP <- calc.indiv.TPs(pop, m.Diag) #calculate new individualised transition probabilities for onset of BC and OC mortality

    m.p <- Probs(m.M[, t], TP) #calculate transition probabilities for 4 states at cycle t (excludes TP for those with invasive BC)
    
    m.M[, t+1] <- samplev(m.p, m.Rand, t) #Sample the next health state and store it in the Matrix m.M numerically
    # 1- no cancer, 2- LG cancer, 3 -HG cancer, 4 - mortality
    
    # Record the characteristics of onset for HG cancer
    m.Diag[, "HG_yr_onset"][m.Diag[, "HG_yr_onset"] >=1] <- m.Diag[, "HG_yr_onset"][m.Diag[, "HG_yr_onset"] >=1] +1 #Update the year of onset if the cancer developed the previous years
    new_HG <-  1*(m.M[, t+1] ==3 & m.M[, t] != 3)
    m.Diag[, "HG_yr_onset"] <- m.Diag[, "HG_yr_onset"] + new_HG
    
    #Mark those individuals who just had BC onset
    m.Diag[, "HG_age_onset"] <- m.Diag[, "HG_age_onset"] + (pop[, "age"] * new_HG)  
    
    # Mark in m.Diag all persons with BC (independently on diagnosis)
    m.Diag[ ,"HG_state"] <- replace(m.Diag[ ,"HG_state"], m.Diag[ ,"HG_yr_onset"] >0, 1)
    
    # Mark age of those individuals who just had LG onset
    new_LG <- m.M[ ,t+1]==2 & m.M[ ,t]==1
    
    # Record the characteristics of onset for LG cancer
    m.Diag[, "LG_state"][which(new_LG)] <- 1
    m.Diag[, "LG_age_onset"] <- m.Diag[, "LG_age_onset"] + (pop[, "age"] * new_LG) 

    # Update m.M_8s matrix with the states for those who have HG cancer
    m.M_8s[, t+1] <- m.M[, t+1]
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.M_8s[, t] ==8, 8) #Replace with BC death those who died with BC before this cycle
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], ((m.M[, t+1]==3 & m.BC.T.to.Stage[ ,"T.onsetToStage2"]==m.Diag[,"HG_yr_onset"])| (m.M[, t+1]==3 & m.M_8s[, t]==5)) & m.Diag[, "HG_yr_diag"] ==0, 5) #Replace for stage 2
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], ((m.M[, t+1]==3 & m.BC.T.to.Stage[ ,"T.onsetToStage3"]==m.Diag[,"HG_yr_onset"])| (m.M[, t+1]==3 & m.M_8s[, t]==6)) & m.Diag[, "HG_yr_diag"] ==0, 6) #Replace for stage 3
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], ((m.M[, t+1]==3 & m.BC.T.to.Stage[ ,"T.onsetToStage4"]==m.Diag[,"HG_yr_onset"])| (m.M[, t+1]==3 & m.M_8s[, t]==7)) & m.Diag[, "HG_yr_diag"] ==0, 7) #Replace for stage 4
    
    # replace with the stage for those who were diagnosed assuming that they don't progress
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "HG_stage_diag"]==2 & m.M[, t+1] != 4, 5) #replace the stage at diagnosis for those who were diagnosed
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "HG_stage_diag"]==3 & m.M[, t+1] != 4, 6) #replace the stage at diagnosis for those who were diagnosed
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], m.Diag[, "HG_stage_diag"]==4 & m.M[, t+1] != 4, 7) #replace the stage at diagnosis for those who were diagnosed

    
    # Update the matrix with the health state numbers with either 0 or 1 depending if it is equal to the sampled state
    m.State[] <- 0
    for(n in 1:n.s_long) {
      m.State[, n] <- replace(m.State[, n], m.M_8s[, t+1] ==n, 1)
    }
    
    #Symptomatic detection
    m.Diag <- f.symptom(m.Diag, m.State, m.Rand, pop, t, m.M) 
    
    # Screening detection
    if(Dipstick_screen ==1){
      scr.Params <- calc.screen.Params(pop, m.Screen, m.State)
    }
  
    # define BC deaths
    TP_BC.mort <- calc.BCmort.TP(pop, m.Diag)
    new_BC1_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.1.mort"]) & m.State[ ,"St1_HG"]>0)
    new_BC2_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.2.mort"]) & m.State[ ,"St2_HG"]>0)
    new_BC3_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.3.mort"]) & m.State[ ,"St3_HG"]>0)
    new_BC4_death <- 1*((m.Rand[ ,"Death_BC", t] < TP_BC.mort[,"TP.BC.4.mort"]) & m.State[ ,"St4_HG"]>0)
    BC_death_all <- rowSums(cbind(new_BC1_death, new_BC2_death, new_BC3_death, new_BC4_death))
    
    # update with the age of BC death
    m.Diag[, "age_BC_death"] <- m.Diag[, "age_BC_death"] + (pop[, "age"] * BC_death_all) # Record the age of death for those with cancer
    
    # Update the mortality for BC (replace with the state 8 in the m.M_8s matrix and state 4 in m.M matrix those who died with BC)
    m.M_8s[, t+1] <- replace(m.M_8s[, t+1], BC_death_all>0, 8) # If died from BC replace to state 8
    m.M[, t+1] <- replace(m.M[, t+1], m.M_8s[, t+1]==8, 4)
    
    # Update the QALYs and costs (not half cycle corrected)
    m.E[, t+1] <- calc.utility(m.State, m.Diag, pop, t) # Assess effects per individual during the cycle t+1 
    m.C[, t+1] <- calc.cost(m.State, m.Diag, m.Cost.treat) # Assess CRC treatment costs per individual during the cycle t+1 (note not half cycle corrected) 
    
    # Gather outcomes for cycle t (total and subgroups)
    m.Out <- aggregate.outcomes(m.M_8s, m.Out, m.M, m.E, m.C, m.Diag, m.Screen, m.State, pop, t)
   
    # Update the age for only those individuals who are alive
    IND <- m.M[, t+1] != 4
    pop[IND, "age"] <- pop[IND,"age"] +1 # update the age
    rand.quit <- m.Rand[ ,"Smoke_quit", t]
    pop[IND, ] <- f.smoke.change(pop[IND, ], rand.quit[IND])  # update smoking status
    pop[IND, ] <- f.risk.calc(pop[IND, ]) #update the risk of BC and p of onset of BC
    
  }#this is a loop for the time point
  

  # Extract the matrices for technical validity
  if(run_mode == "Testing") { # Create a matrix of transitions across states
    TS_8 <- paste(m.M_8s, cbind(m.M_8s[, -1], NA), sep = "->") # transitions from one state to the other     
    TS_8 <- matrix(TS_8, nrow = n.i)
    TS_4 <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other     
    TS_4 <- matrix(TS_4, nrow = n.i)
    rownames(TS_8) <- rownames(TS_4) <-paste("Ind",   1:n.i, sep = " ")   # name the rows      
    colnames(TS_8) <- colnames(TS_4) <-paste("Cycle", 0:n.t, sep = " ")   # name the columns    
  
    TR_8 <- t(apply(m.M_8s, 2, function(x) table(factor(x, levels = v.n_long, ordered = TRUE))))     
    TR_8 <- TR_8 / n.i    # create a distribution trace
    TR_4 <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))     
    TR_4 <- TR_4 / n.i    # create a distribution trace
    
    rownames(TR_8) <- rownames(TR_4) <-paste("Cycle", 0:n.t, sep = " ")  # name the rows
    colnames(TR_8) <- v.n_long     # name the columns    
    colnames(TR_4) <-v.n    # name the columns    
    
  } else {   
    TS_8 <- TS_4 <- TR_8 <- TR_4 <- NULL   
  } 
  
  # Extract the matrices by sex for calibration
  
  if (run_mode == "Calibration") { #  add additional two matrices to report TR matrix by gender
    names_Diag <- c(colnames(m.Diag), "sex", "age")
    m.Diag <- cbind(m.Diag, pop[ ,"sex"], pop[ ,"age_0"])
    colnames(m.Diag) <-names_Diag
    
    a.Temp <- array(data = NA, dim = c(n.i, n.t+1, n.s_long))
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
  
  results <- switch(run_mode, 
                    Testing = list(TR_8 = TR_8, TR_4 = TR_4,m.Diag = m.Diag, m.State=m.State, m.Out=m.Out),
                    Calibration = list(TR_f = TR_f, TR_m=TR_m, m.Diag = m.Diag),
                    Deterministic = m.Out,
                    PSA= m.Out
                    )
  
  
  return(results)  # return the results 
  }  