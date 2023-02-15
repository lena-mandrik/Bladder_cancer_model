# This is the code to run calibration either with a random approach or a Bayesian way
# Date: 23/11/22

rm(list = ls(all = TRUE)) # Remove everything in the workspace 

#Load required libraries
library("survival")
library("survminer")
library("dplyr")
library("vroom")
library("purrr")
library("foreach")
library("doParallel")
library("doRNG")
library(lhs) # for Random search using Latin-Hypercube Sampling

###Set random seed
set.seed(10)


## Load all data needed for the model that is not in the loop
##################################################################################

###Set up the Global Parameters
run_mode <- "Calibration" # Available modes include "Testing" (returns all matrices), "Calibration_rand" (m.Diag and TR), "Calibration_Bayes" (m.Diag and TR),"Deterministic" (m.Out only), "PSA" (m.Out only)
calibration_type = "Bayes" #also can be Random
cohort <- 1 # 1 = all individuals start model at same age (cohort), 0 = individuals start in model at true (HSE) age
cohort_age <- 30 #select starting age of cohort (hash out or set to anything if not using cohort)
n.loops <- 3 # The number of model loops/PSA loops to run 
cl <- 1  # The cycle length (years) 
n.t   <- if(cohort==1){100-cohort_age}else{70}  # The number of cycles to run 
d.c <- 0.035 # The discount rate for costs
d.e <- 0.035 # The discount rate for effects
N_sets <- if(calibration_type =="Bayes" ){n.loops}else{1} #Number of parameter sets required (for PSA). Minimum value = 1 (mean parameter values); maximum value = 1651 (number of calibrated correlated param sets)

###Load up all the functions and all the data for use in the model
source("Load_all_files.R")


###Load up all the functions and all the data for use in the model
###############################################################################
# Load the calibration inputs
# Read all functions from all scripts within the R calibration folder

#sapply(list.files("R calibration/R",full.names=T), source)
source("R calibration/R/func_calibration.R")

# Calculate the size of the cohorts 
Cohort_m <- sum(population[,"sex"]==1)
Cohort_f <- sum(population[,"sex"]==0)

# Extract the calibration targets and calculate the confidence interval
# The incidence is already converted per 1 person alive
Targets <- read.table("R calibration/Calibration_targets.txt", header = TRUE, row.names=1) # Load targets by 5-years
Targets <- f.targets.per.alive.long(Targets, start_age=30, end_age =95, n_by=5) #Linearly extrapolate between the 5-year bands

n_targets <- ncol(Targets) -1 # Exclude Age column
v.target_names <- colnames(Targets)[-1] # Exclude Age column

# USE CI with beta distribution
CI_targets <- read.table("R calibration/CI_targets.txt", header = TRUE) # Beta distribution for the CI # f.CI.calc(Targets)
CI_targets <- f.targets.per.alive.long(CI_targets, start_age=30, end_age =95, n_by=5) #Linearly extrapolate between the 5-year bands


# Calculate SE from the CI
SE_males <- (CI_targets[ ,14:19]- CI_targets[ ,2:7])/3.92
SE_females <- (CI_targets[ ,20:25] - CI_targets[ ,8:13])/3.92
SE <- cbind(SE_males, SE_females)
colnames(SE) <- paste(colnames(Targets)[-1],"_SE", sep="")

# Visualize the target
f.plot.target(Targets, CI_targets)

# Define parameters to calibrate
Calibr_parameters <- Params[c("P.onset", "P.onset_low.risk", "P.onset_age", "RR.onset_sex", "P.sympt.diag_LGBC", "P.sympt.diag_A_HGBC",
                              "P.sympt.diag_B_HGBC", "P.sympt.diag_Age80_HGBC",
                               "shape.t.StI.StII", "shape.t.StII.StIII", "shape.t.StIII.StIV", "P.LGtoHGBC") ,]
v.param_names <- rownames(Calibr_parameters) # number of parameters to calibrate 
n_params <- length(v.param_names)

n_samples <-2 # number of calibration runs

# Goodness-of-fit outcomes
m.GOF <-  matrix(nrow=n_samples, ncol = ncol(Targets)-1)
colnames(m.GOF) <- paste0(c(v.target_names), "_fit")

if(run_mode == "Calibration_rand"){
  source("R calibration//R//Random//f.random.R")
  print(output)
}