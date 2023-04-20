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
run_mode <- "PSA" # Available modes include "Testing" (returns all matrices), "Calibration" (m.Diag and TR), "Deterministic" (m.Out only), "PSA" (m.Out only)
cohort <- 1 # 1 = all individuals start model at same age (cohort), 0 = individuals start in model at true (HSE) age
cohort_age <- 30 #select starting age of cohort (hash out or set to anything if not using cohort)
n.loops <- 200 # The number of model loops/PSA loops to run 
cl <- 1  # The cycle length (years) 
n.t   <- if(cohort==1){100-cohort_age}else{70}  # The number of cycles to run 
d.c <- 0.035 # The discount rate for costs
d.e <- 0.035 # The discount rate for effects
N_sets <- if(run_mode =="PSA"){n.loops}else{1} #Number of parameter sets required (for PSA). Minimum value = 1 (mean parameter values); maximum value = 1651 (number of calibrated correlated param sets)

###Load up all the functions and all the data for use in the model
source("Load_all_files.R")

###Load up all the functions and all the data for use in the model
###Set up results collection
results_no_screen <- results_screen <-list()

# create a list to collect PSA outcomes if the model runs in the PSA mode
if(run_mode == "PSA"){
  PSA_results_no_screen <- matrix(0, nrow = 8 * length(out_names), ncol = n.loops)
  rownames(PSA_results_no_screen) <- rep(out_names, 8)
  PSA_results_screen <- PSA_results_no_screen
}
### run the model:
system.time(for(iter in 1:n.loops) {
  
#Select appropriate parameter set to use
p.set <- ifelse(run_mode =="PSA", iter, 1) 
#p.set=1
#Set up the model parameters according to current parameter set
f.set_parameters(p.set)

# Allocate the time to stage at diagnosis for each person in HSE
m.BC.T.to.Stage <- f.stage(Mean.t.StI.StII, shape.t.StI.StII, Mean.t.StII.StIII, shape.t.StII.StIII, 
                           Mean.t.StIII.StIV, shape.t.StIII.StIV, n.i)

#calculates relative and absolute risks of cancer onset
pop <- f.risk.calc(population) 

#Set up random number array for each individual
m.Rand <- f.generate_random()

Dipstick_screen =0 #Set whether the screening with dipstick happens, 0 - no, 1- yes
Dipstick_age =0 #Set the age of the dipstick if screening happens, set to zero if no screening
results_no_screen[[iter]] = Simulate_NHD(n.i, n.t, pop)

if(run_mode == "PSA"){
  PSA_results_no_screen[, iter] <- rowSums(results_no_screen[[iter]])
  write.table(PSA_results_no_screen,"PSA_results_no_screen.txt")
}

#Dipstick_screen =0 #Set whether the screening with dipstick happens, 0 - no, 1- yes
#Dipstick_age =100 #Set the age of the dipstick if screening happens
#results_screen[[iter]] = Simulate_NHD(n.i, n.t, pop)

})

############################################
#' Analyse the results depending whether the runds were deterministic or probabilistic 
#' Analyse the results
if(run_mode == "Deterministic"){
  
  results_no_screen_cum <- process_DA_results(results_no_screen, "results_NoScreen.txt")
  
}

if(run_mode == "PSA"){
  
  results_no_screen_PSA <- process_PSA_results(results_no_screen, "results_no_screen_PSA.txt")
}