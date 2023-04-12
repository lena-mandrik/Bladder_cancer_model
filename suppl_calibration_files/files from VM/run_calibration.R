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
set.seed(333)


## Load all data needed for the model that is not in the loop
##################################################################################

###Set up the Global Parameters
run_mode <- "Calibration" # Available modes include "Testing" (returns all matrices), "Calibration" (m.Diag and TR), "Deterministic" (m.Out only), "PSA" (m.Out only)
cohort <- 1 # 1 = all individuals start model at same age (cohort), 0 = individuals start in model at true (HSE) age
cohort_age <- 30 #select starting age of cohort (hash out or set to anything if not using cohort)
n.loops <- 100 # The number of model loops/PSA loops to run 
cl <- 1  # The cycle length (years) 
n.t   <- if(cohort==1){100-cohort_age}else{70}  # The number of cycles to run 
d.c <- 0.035 # The discount rate for costs
d.e <- 0.035 # The discount rate for effects
N_sets <- if(run_mode =="PSA"){n.loops}else{1} #Number of parameter sets required (for PSA). Minimum value = 1 (mean parameter values); maximum value = 1651 (number of calibrated correlated param sets)

###Load up all the functions and all the data for use in the model
source("Load_all_files.R")


###Load up all the functions and all the data for use in the model
###############################################################################
# Load the calibration inputs

# Read all functions from all scripts within the R calibration folder

sapply(list.files("R calibration/R",full.names=T), source)

# Select the type of the calibration
calibr_type <- "Random" # The options will include random and Bayesian ways

# Calculate the size of the cohorts 
Cohort_m <- sum(population[,"sex"]==1)
Cohort_f <- sum(population[,"sex"]==0)

# Extract the calibration targets and calculate the confidence interval
# The incidence is already converted per 1 person alive
Targets <- read.table("R calibration/Calibration_targets.txt", header = TRUE, row.names=1) # Load targets by 5-years
Targets <- f.targets.per.alive.long(Targets, start_age=30, end_age =95, n_by=5) #Linearly extrapolate between the 5-year bands

n_targets <- ncol(Targets) -1 # Exclude Age column
v.target_names <- colnames(Targets)[-1] # Exclude Age column
CI_targets <- f.CI.calc(Targets)

# Calculate SE from the CI
SE_males <- (CI_targets[ ,7:12] - CI_targets[ ,1:6])/3.92
SE_females <- (CI_targets[ ,19:24] - CI_targets[ ,13:18])/3.92
SE <- cbind(SE_males, SE_females)
colnames(SE) <- paste(colnames(Targets)[-1],"_SE", sep="")

# Visualize the target
f.plot.target(Targets, CI_targets)

# Define parameters to calibrate
Calibr_parameters <- Params[c("P.onset", "P.onset_low.risk", "P.onset_age", "RR.onset_sex", "P.sympt.diag_LGBC", "P.sympt.diag_A_HGBC",
                              "P.sympt.diag_B_HGBC", "P.sympt.diag_Age80_HGBC", "C.age.80.undiag.mort",
                               "shape.t.StI.StII", "shape.t.StII.StIII", "shape.t.StIII.StIV") ,]
v.param_names <- rownames(Calibr_parameters) # number of parameters to calibrate 
n_params <- length(v.param_names)

n_samples <-5000 # number of calibration runs

time_0 <- Sys.time() # Add in the end the run time calculation

# Search dimensions. Use for random search algorithm only
lower_boud <- c(0.000001, 0.1, 1.001, 1.001, 0.005, 0.001, 1.01, 0.85, -0.3, rep(1.005,3)) #for symtomatic presentation requires other limits - A <15%, B =5-40%, C=10-60%, D-30-100%
upper_boud <- c(0.001, 0.6, 1.2, 1.9, 0.1, 0.2, 1.9, 1, -0.01, rep(1.5,3))

# Sample using Latin Hypercube
sample_LH <- randomLHS(n_samples, n_params)

# Re-scale to min/max
m.sample.params <- matrix(nrow=n_samples, ncol = n_params)
colnames(m.sample.params) <- v.param_names

for(i in 1:n_params){
  m.sample.params[,i] <- qunif(sample_LH[,i],
                               min = lower_boud[i],
                               max = upper_boud[i])
}
write.csv(m.sample.params, file="R calibration\\ParametersInput.csv")

# Goodness-of-fit outcomes
m.GOF <-  matrix(nrow=n_samples, ncol = ncol(Targets)-1)
colnames(m.GOF) <- paste0(c(v.target_names), "_fit")


# Run the calibration loops

for(run in 1:n_samples){
  
  Calibr_parameters <- as.matrix(m.sample.params[run,], ncol = 1) # get the inputs for calibrating parameters
  
  source("R calibration\\Main_model_calibration.R") # run the model
  
  output <- f.calibr.output(results_no_screen) # extract the outputs
  
  Predict <- cbind(output$rate_outcomes_m, output$rate_outcomes_f)
  
  # log likelihood instead of sum of squared errors
  
  m.GOF <- f.GOF.calc(run, m.GOF, Targets, SE, Predict)
  

  if(run%%10==0){
    cat('\r', paste(round(run/n_samples*100), "% done", sep = " "))
    write.csv(m.GOF, file="R calibration\\Outputs\\m.GOF.csv")}
} 

