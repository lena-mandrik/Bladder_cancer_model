# This is the code to run calibration either with a random approach or a Bayesian way

# Select the type of the calibration

calibr_type <- "Random" # The options will include random and Bayesian ways


# Calculate the size of the cohorts 
Cohort_m <- sum(population[,"sex"]==1)
Cohort_f <- sum(population[,"sex"]==0)


# Extract the calibration targets
Targets <- as.matrix(read.table("Data/Calibration_targets.txt", header = TRUE, row.names=1)) # Load targets by 5-years

Targets <- f.targets.per.alive.long(Targets, start_age=30, end_age =95, n_by=5) #Linearly extrapolate between the 5-year bands

n_targets <- ncol(Targets) -1 # Exclude Age column
v.target_names <- colnames(Targets)[-1] # Exclude Age column

CI_targets <- f.CI.calc(Targets)

  

  
# Visualize the target
targets_pt = ggplot(data = lst_targets$Surv,
                    aes(x = time,
                        y = value)) +
  geom_errorbar(aes(ymin = lb, ymax = ub)) +
  geom_point() +
  theme(
    panel.border = element_rect(fill = NA, color = 'black')
  ) + 
  labs(title = "Calibration target",
       x = "Time",
       y = "Proportion survived")
targets_pt

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