# Lena Mandrik 12.07.2022
# This is the main script
###Run this script to run whole model
###Remove everything in workspace 

rm(list = ls(all = TRUE))

#Load required libraries


###Set random seed
set.seed(10)


###Set up the Global Parameters
run_mode <- "Deterministic" 
cohort <- 1 # 1 = all individuals start model at same age (cohort), 0 = individuals start in model at true (HSE) age
cohort_age <- 30 #select starting age of cohort (hash out or set to anything if not using cohort)

###Load up all the functions and all the data for use in the model

## Here the loops will be

#Select appropriate parameter set to use
p.set <- ifelse(run_mode =="PSA", i, 1) 