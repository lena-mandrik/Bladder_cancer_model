### Function to run the simulation ###

#Simulate_NHD <- function(n.i, n.t, pop, FITscreen, FITage, FSscreen, FSage, FITthresh, model, risk, sex) { 
  # Arguments:     
  # n.i:     number of individuals; Global parameters, Master script
  # n.t:     total number of cycles to run the model; Global parameters, Master script   
  # pop:     matrix of individual level population attributes
  # FITscreen:  whether or not FIT screening is used in this model run
  # FITage: age at which FIT screening starts (if selected)
  # FSscreen:  whether or not FS screening is used in this model run
  # FSage: age at which one of FS screen is given
  # model: model used to calculate individual CRC risk
  # risk: absolute risk at which FIT screening starts if using risk model
  # sex: use sex-specific risks to calculate age at first screen
  
  # Makes use of:   
  # Probs:  function for the estimation of transition probabilities (script "Functions")
  # Costs:   function for the estimation of cost state values (script "Functions")
  # Effs:    function for the estimation of state specific health outcomes (QALYs) (script "Functions")   
  
  # create matrices capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.C <- m.E <-  matrix(nrow = n.i, ncol = n.t + 1)
  
  m.M[, 1] <- rep(1, n.i)   # indicate the initial health state
  m.C[, 1] <- 0 # estimate costs per individual for the initial health state 
  m.E[, 1] <- pop[, "EQ5D"] - (Utility_age * (pop[, "age"] - pop[, "age_0"])) # estimate QALY per individual for the initial health state and starting age
  m.E[, 1] <- replace(m.E[, 1], m.E[, 1] >1, 1) # if EQ-5D over 1 reset to 1
  m.E[, 1] <- replace(m.E[, 1], m.E[, 1] <=(-0.594), -0.594) # if EQ-5D under -0.594 reset to -0.594
  
  # create another matrix capturing the current health state for each individual
  m.State <- matrix(0, nrow = n.i, ncol = n.s)
  colnames(m.State) <- states
  for(n in 1:n.s) {
    m.State[, n] <- replace(m.State[, n], m.M[, 1] ==n, 1)
  }
  
  #Create another matrix for current diagnostic information
  m.Diag <- matrix(0, nrow = n.i, ncol = 8)
  colnames(m.Diag) <- c("BC_diag", "sympt_diag", "screen_diag", "new_diag",
                        "age_diag", "stage_diag", "yr_onset", "yr_diag")
  
  #Create an array to gather screening and surveillance information for each cycle
  m.Screen <- array(data = 0, dim = c(n.i, n.t + 1, length(screen_names)))
  dimnames(m.Screen)[[3]] <- screen_names # note that this is defined in set_params
  
  #Create a matrix for gathering aggregate outcomes
  m.Out <- matrix(0, nrow = length(out_names), ncol = n.t + 1)
  rownames(m.Out) <- out_names
  colnames(m.Out) <- c(0:n.t)
  
  
  # loop to run the model over time
  for(t in 1:n.t) {
    
    #Natural History
    TP <- calc.indiv.TPs(pop, m.Diag, m.State) #calculate new individualised transition probabilities for each health state
    
    m.p <- Probs(m.M[, t], TP) #calculate transition probabilities at cycle t
    
    m.M[, t+1] <- samplev(m.p, m.Rand, t) #Sample the next health state and store it in the Matrix m.M numerically
    
    # Update the matrix with the health state numbers with either 0 or 1 depending if it is equal to the sampled state
    m.State[] <- 0
    for(n in 1:n.s) {
      m.State[, n] <- replace(m.State[, n], m.M[, t+1] ==n, 1)
    }
    
    #Symptomatic detection
    m.Diag <- f.symptom(m.Diag, m.State, m.Rand, pop, t)
    