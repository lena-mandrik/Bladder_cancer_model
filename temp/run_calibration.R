# This is the code to run calibration either with a random approach or a Bayesian way
# Date: 23/11/22


# Read all functions from all scripts within the R calibration folder

sapply(list.files("R calibration",full.names=T), source)

# Select the type of the calibration

calibr_type <- "Random" # The options will include random and Bayesian ways


# Calculate the size of the cohorts 
Cohort_m <- sum(population[,"sex"]==1)
Cohort_f <- sum(population[,"sex"]==0)


# Extract the calibration targets and calculate the confidence interval
# The incidence is already converted per 1 person alive
Targets <- read.table("Data/Calibration_targets.txt", header = TRUE, row.names=1) # Load targets by 5-years
Targets <- f.targets.per.alive.long(Targets, start_age=30, end_age =95, n_by=5) #Linearly extrapolate between the 5-year bands

n_targets <- ncol(Targets) -1 # Exclude Age column
v.target_names <- colnames(Targets)[-1] # Exclude Age column
CI_targets <- f.CI.calc(Targets)


# Visualize the target
f.plot.target(Targets, CI_targets)


# Define parameters to calibrate
Calibr_parameters <- Params[c("P.onset", "P.onset_age", "P.onset_sex", "P.sympt.diag_LGBC", "P.sympt.diag_A_HGBC",
                              "P.sympt.diag_B_HGBC", "P.sympt.diag_Age80_HGBC", "C.age.80.undiag.mort",
                              "RR.All.Death.no_smoke", "shape.t.StI.StII", "shape.t.StII.StIII", "shape.t.StIII.StIV") ,]
v.param_names <- rownames(Calibr_parameters) # number of parameters to calibrate 
n_params <- length(v.param_names)

n_samples <-8 #number of calibration runs

# Extracting from the lists 
list_results <- results_no_screen
list <- f.calibr.output(list_results)