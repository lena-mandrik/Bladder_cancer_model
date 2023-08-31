# Lena Mandrik 12.07.2022
# This is the main script
###Run this script to run whole model


# Random numbers:
optsN <- list(123, normal.kind = "Ahrens")

# Linux or Windows code to forking or parallel execution:
if(f.get_os() == "windows") {
  n_cores <- (detectCores()-2)
  cluster = makeCluster(n_cores, type = "PSOCK", outfile = "")  # Windows - WM
  registerDoParallel(cluster) # Windows- WM
  #registerDoRNG(1234, once = FALSE) # Windows- WM
} else {
  #doMC::registerDoMC(10) # Linux - WM. Random number issues
  cluster = makeCluster(10, type = "FORK", outfile = "")  # Linux - WM
  registerDoParallel(cluster) # Linux- WM
  #registerDoRNG(1234, once = FALSE) # Linux- WM
  # Progress bar:
  pb = txtProgressBar(min = 1, max = n.loops, style = 3) # Linux- WM
}



### run the model:
results_no_screen = foreach::foreach(iterator = 1:n.loops, .options.RNG = optsN) %dorng% {
  gc()
  
  #Progress bar:
  if(f.get_os() != "windows") {
    setTxtProgressBar(pb, iterator) 
  }
  
#Select appropriate parameter set to use
p.set <- ifelse(calibration_type =="Bayes", iterator, 1)


#Set up the model parameters according to current parameter set
f.set_parameters(p.set)

P.Recurrence.LR <-0 #set recurrence to zero

  P.onset <<- Calibr_parameters["P.onset", 1]
  P.onset_low.risk <<- Calibr_parameters["P.onset_low.risk", 1]
  P.onset_age <<- Calibr_parameters["P.onset_age", 1]
  RR.onset_sex <<- Calibr_parameters["RR.onset_sex", 1]
  
  P.sympt.diag_LGBC <<- Calibr_parameters["P.sympt.diag_LGBC", 1]
  P.sympt.diag_St1 <<- Calibr_parameters["P.sympt.diag_St1", 1]
  P.sympt.diag_St2 <<- Calibr_parameters["P.sympt.diag_St2", 1]
  P.sympt.diag_St3 <<- Calibr_parameters["P.sympt.diag_St3", 1]
  P.sympt.diag_St4 <<- Calibr_parameters["P.sympt.diag_St4", 1]
  
  
  P.sympt.diag_Age <<- Calibr_parameters["P.sympt.diag_Age", 1]

  shape.t.StI.StII <<- Calibr_parameters["shape.t.StI.StII", 1]
  shape.t.StII.StIII <<- Calibr_parameters["shape.t.StII.StIII", 1]
  shape.t.StIII.StIV <<- Calibr_parameters["shape.t.StIII.StIV", 1]
  P.LGtoHGBC <<- Calibr_parameters["P.LGtoHGBC", 1]
  P.ungiag.dead <<- Calibr_parameters["P.ungiag.dead", 1]
  
  
  ######### increase mean time
  Mean.t.StI.StII <<- Calibr_parameters["Mean.t.StI.StII", 1] #5.93
  Mean.t.StII.StIII <<- Calibr_parameters["Mean.t.StII.StIII", 1] #4.42
  Mean.t.StIII.StIV <<- Calibr_parameters["Mean.t.StIII.StIV", 1] #2.7
  
# Allocate the time to stage at diagnosis for each person in HSE
m.BC.T.to.Stage <- f.stage(Mean.t.StI.StII, shape.t.StI.StII, Mean.t.StII.StIII, shape.t.StII.StIII, 
                           Mean.t.StIII.StIV, shape.t.StIII.StIV, n.i)

#calculates relative and absolute risks of cancer onset
pop <- set_population(population) 

#Set up random number array for each individual
m.Rand <- f.generate_random()

# Update the smoking status for the previous smokers with the current smokers
pop <- f.smoke.initial(pop, m.Rand)

Simulate_NHD(nsample, n.t, pop, m.BC.T.to.Stage)

}

if(f.get_os() == "windows") {
  stopCluster(cluster)
} else {
  close(pb)
  stopCluster(cluster)
}
