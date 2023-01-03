# Lena Mandrik 12.07.2022
# This is the main script
### Run this script to run whole model
### Remove everything in workspace 

rm(list = ls(all = TRUE))

#Load required libraries
library("survival")
library("survminer")
library("dplyr")
library("vroom")
library("purrr")
library("foreach")
library("doParallel")
library("doRNG")

###Set random seed
set.seed(10)


###Set up the Global Parameters
run_mode <- "Testing" # Available modes include "Testing" (returns all matrices), "Calibration" (m.Diag and TR), "Deterministic" (m.Out only), "PSA" (m.Out only)
cohort <- 1 # 1 = all individuals start model at same age (cohort), 0 = individuals start in model at true (HSE) age
cohort_age <- 30 #select starting age of cohort (hash out or set to anything if not using cohort)
n.loops <- 1 # The number of model loops/PSA loops to run 
cl <- 1  # The cycle length (years) 
n.t   <- if(cohort==1){100-cohort_age}else{70}  # The number of cycles to run 
d.c <- 0.035 # The discount rate for costs
d.e <- 0.035 # The discount rate for effects
N_sets <- if(run_mode =="PSA"){n.loops}else{1} #Number of parameter sets required (for PSA). Minimum value = 1 (mean parameter values); maximum value = 1651 (number of calibrated correlated param sets)

###Load up all the functions and all the data for use in the model
source("Load_all_files.R")


###Load up all the functions and all the data for use in the model
###Set up results collection
results_no_screen <- list()

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
p.set <- ifelse(run_mode =="PSA", i, 1) 

#Set up the model parameters according to current parameter set
f.set_parameters(p.set)

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

