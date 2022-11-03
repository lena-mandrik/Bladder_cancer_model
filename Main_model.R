# Lena Mandrik 12.07.2022
# This is the main script
###Run this script to run whole model
###Remove everything in workspace 

rm(list = ls(all = TRUE))

#Load required libraries
library("survival")
library("survminer")
library("dplyr")
library("vroom")
library("purrr")

###Set random seed
set.seed(10)


###Set up the Global Parameters
run_mode <- "Deterministic" # Available modes include "Testing" (returns all matrices), "Calibration" (m.Diag only), "Deterministic" (m.Out only), "PSA" (m.Out only)
cohort <- 1 # 1 = all individuals start model at same age (cohort), 0 = individuals start in model at true (HSE) age
cohort_age <- 50 #select starting age of cohort (hash out or set to anything if not using cohort)
n.loops <- 5 # The number of model loops/PSA loops to run 
cl <- 1  # The cycle length (years) 
n.t   <- if(cohort==1){100-cohort_age}else{70}  # The number of cycles to run 
d.c <- 0.035 # The discount rate for costs
d.e <- 0.035 # The discount rate for effects
N_sets <- if(run_mode =="PSA"){n.loops}else{1} #Number of parameter sets required (for PSA). Minimum value = 1 (mean parameter values); maximum value = 1651 (number of calibrated correlated param sets)


#score_age <- 40  #(Default = 40). Age at which known family history assigned and CRC risk score is calculated. 
#CRC_death <- 1 #(Default = 1). Multiplier to adjust annual probability of CRC death. Globally applied to all ages, stages, sexes and times since diagnosis.
#CRC_inc <- 1 #(Default = 1). Multiplier to adjust individual CRC risk and thereby incidence. Globally applied to all individuals.

###Load up all the functions and all the data for use in the model
source("Load_all_files.R")

#Set up results collection for all simulations named below
results_NoScreen <- results_screen <- list()

#Set up PSA results collection
if(run_mode == "PSA"){
  PSA_results_NoScreen <- matrix(0, nrow = 8 * length(out_names), ncol = n.loops)
  rownames(PSA_results_NoScreen) <- rep(out_names, 8)
  PSA_results_screen <- PSA_results_NoScreen
}

###Load up all the functions and all the data for use in the model
## Here the loops will be
i=1

#Select appropriate parameter set to use
p.set <- ifelse(run_mode =="PSA", i, 1) 

#Set up the model parameters according to current parameter set
set_parameters(p.set)

# Allocate the time to stage at diagnosis for each person in HSE
m.BC.T.to.Stage <- f.stage.assign(m.BC.T.to.Stage) ##!! Needs to be replaced later to avoid loops (apply or map)

# Set the model baseline population risk
pop <- f.risk.calc(population)

# Something else 