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
p.set <- 1

#Set up the model parameters according to current parameter set
f.set_parameters(p.set)

if(run_mode == "Calibration_rand"){ #replace the calibrated parameters in the calibration mode only with the sampled deterministic values in random calibration
  P.onset_age <- Calibr_parameters["P.onset_age", 1]
  P.onset_low.risk <- Calibr_parameters["P.onset_low.risk", 1]
  P.onset_sex <- Calibr_parameters["RR.onset_sex", 1]
  P.sympt.diag_LGBC <- Calibr_parameters["P.sympt.diag_LGBC", 1]
  P.sympt.diag_A_HGBC <- Calibr_parameters["P.sympt.diag_A_HGBC", 1]
  P.sympt.diag_B_HGBC <- Calibr_parameters["P.sympt.diag_B_HGBC", 1]
  P.sympt.diag_Age80_HGBC <- Calibr_parameters["P.sympt.diag_Age80_HGBC", 1]
  C.age.80.undiag.mort<- Calibr_parameters["C.age.80.undiag.mort", 1]
  shape.t.StI.StII <- Calibr_parameters["shape.t.StI.StII", 1]
  shape.t.StII.StIII <- Calibr_parameters["shape.t.StII.StIII", 1]
  shape.t.StIII.StIV <- Calibr_parameters["shape.t.StIII.StIV", 1]
}



if(run_mode =="Calibration_Bayes"){ #replace the calibrated parameters in the calibration mode only with values sampled from the distribution
  P.onset_age <- Calibr_parameters["P.onset_age", 1]
  P.onset_low.risk <- Calibr_parameters["P.onset_low.risk", 1]
  P.onset_sex <- Calibr_parameters["RR.onset_sex", 1]
  P.sympt.diag_LGBC <- Calibr_parameters["P.sympt.diag_LGBC", 1]
  P.sympt.diag_A_HGBC <- Calibr_parameters["P.sympt.diag_A_HGBC", 1]
  P.sympt.diag_B_HGBC <- Calibr_parameters["P.sympt.diag_B_HGBC", 1]
  P.sympt.diag_Age80_HGBC <- Calibr_parameters["P.sympt.diag_Age80_HGBC", 1]
  C.age.80.undiag.mort<- Calibr_parameters["C.age.80.undiag.mort", 1]
  shape.t.StI.StII <- Calibr_parameters["shape.t.StI.StII", 1]
  shape.t.StII.StIII <- Calibr_parameters["shape.t.StII.StIII", 1]
  shape.t.StIII.StIV <- Calibr_parameters["shape.t.StIII.StIV", 1]
  
 # x <- runif(50, min=Calibr_parameters[1,]*0.9, max=Calibr_parameters[1,]*1.1)
}


# Allocate the time to stage at diagnosis for each person in HSE
m.BC.T.to.Stage <- f.stage(Mean.t.StI.StII, shape.t.StI.StII, Mean.t.StII.StIII, shape.t.StII.StIII, 
                           Mean.t.StIII.StIV, shape.t.StIII.StIV, n.i)

#calculates relative and absolute risks of cancer onset
pop <- f.risk.calc(population) 

#Set up random number array for each individual
m.Rand <- generate_random()

Simulate_NHD(n.i, n.t, pop)

}

if(f.get_os() == "windows") {
  stopCluster(cluster)
} else {
  close(pb)
  stopCluster(cluster)
}
