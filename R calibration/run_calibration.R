# This is the code to run calibration either with a random approach or a Bayesian way
# Date: 23/11/22

rm(list = ls(all = TRUE)) # Remove everything in the workspace 

#Load required libraries
library("dplyr")
library("ggpubr")
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
set.seed(137)


## Load all data needed for the model that is not in the loop
##################################################################################

###Set up the Global Parameters
run_mode <- "Calibration" # Available modes include "Testing" (returns all matrices), "Calibration" (m.Diag and TR),"Deterministic" (m.Out only), "PSA" (m.Out only)
calibration_type = "Bayes"
cohort <- 1 # 1 = all individuals start model at same age (cohort), 0 = individuals start in model at true (HSE) age
cohort_age <- 30 #select starting age of cohort (hash out or set to anything if not using cohort)
n.loops <- 200 # The number of model loops/PSA loops to run 
cl <- 1  # The cycle length (years) 
n.t   <- if(cohort==1){100-cohort_age}else{70}  # The number of cycles to run 
d.c <- 0.035 # The discount rate for costs
d.e <- 0.035 # The discount rate for effects
N_sets <- if(run_mode =="Calibration" ){n.loops}else{1} #Number of parameter sets required (for PSA). Minimum value = 1 (mean parameter values); maximum value = 1651 (number of calibrated correlated param sets)
P.Recurrence.LR <-0 #set up recurrent cancers to zero as they are not reflected in the statistical data
DS_screen =0 #Set whether the screening with dipstick happens, 0 - no, 1- yes
DS_age = 100 #Set the age of the dipstick if screening happens, set to any if no screening

# Set whether the population needs to be all population, current smokers, or both current and past smokers to be sampled before the model start
char_pop <- "all" # c.smoke- current smokers, all.smoke- current and past smokers, all -the whole population
screen_elig <- "all" # c.smoke- current smokers, all.smoke- current and past smokers, all -the whole population
sex <- "all"

disease ="bladder"

###Load up all the functions and all the data for use in the model
source("Load_all_files.R")
nsample <- if(cohort==1){n.i} else{5000} #define a sample size for population to run in each loop 


###Load up all the functions and all the data for use in the model
###############################################################################
# Load the calibration inputs
# Read all functions from all scripts within the R calibration folder

source("R calibration/R/func_calibration.R")

# Calculate the size of the cohorts 
Cohort_m <- sum(population[,"sex"]==1)
Cohort_f <- sum(population[,"sex"]==0)

# Extract the calibration targets and calculate the confidence interval
# The incidence is already converted per 1 person alive
# Try with counts and rates

target_stat = "rates" # also rates

  if(disease =="bladder"){
  #Targets <- as.matrix(read.table("R calibration/bladder_targets.txt", header =FALSE, sep="\t"))# Load targets by 5-years
  Targets <- as.matrix(read.table("R calibration/Targets/Calibration_targets_bladder_counts.txt", header =T))# Load targets by 5-years
  
  CI_targets <- read.table("R calibration/Targets/CI_targets_bladder_counts.txt", header = TRUE) # Beta distribution for the CI # f.CI.calc(Targets)
  }
  #Targets <- f.targets.per.alive.long(Targets, start_age=30, end_age =90, n_by=5) #Linearly extrapolate between the 5-year bands
  #CI_targets <- f.targets.per.alive.long(CI_targets, start_age=30, end_age =90, n_by=5) #Linearly extrapolate between the 5-year bands
  Targets <- Targets[, -which(colnames(Targets)=="Age_lb")]


 # Load pop by 5-years
England.pop <- as.matrix(read.table("Data/Pop.alive.by.age.txt", header = TRUE)) # Load targets by 5-years

  #extrapolate from 5-year data
#England.pop <- f.targets.per.alive.long(England.pop, start_age=30, end_age =95, n_by=5) #Linearly extrapolate between the 5-year bands

#Get CI for the targets
#CI_targets <- read.table("R calibration/Targets/CI_targets.txt", header = TRUE) # Beta distribution for the CI # f.CI.calc(Targets)
England.pop <-England.pop[1:nrow(CI_targets),]

Targets[,1:7] <-Targets[,1:7]/(England.pop[,"Males"]*1000)
Targets[,8:14] <-Targets[,8:14]/(England.pop[,"Females"]*1000)


CI_targets[ ,c(1:7,15:21)] <-CI_targets[ ,c(1:7,15:21)]/(England.pop[,"Males"]*1000)
CI_targets[ ,c(8:14,22:28)] <-CI_targets[ ,c(8:14,22:28)]/(England.pop[,"Females"]*1000)

#p_lower = Targets -1.96*sqrt(Targets)
#colnames(p_lower) <- paste(colnames(Targets),"_lb", sep="")

#p_higher = Targets +1.96*sqrt(Targets)
#colnames(p_higher) <- paste(colnames(Targets),"_ub", sep="")

#CI_targets <- cbind(p_lower, p_higher)
#CI_targets[CI_targets <0] <-0

# Load pop by 5-years
#England.pop_old <- colSums(England.pop[13:14,])
#England.pop <- rbind(England.pop[1:12,],England.pop_old) 

#n_targets <- ncol(Targets) 
#v.target_names <- colnames(Targets) 


  # USE CI with beta distribution
  # Calculate SE from the CI
  SE_males <- (CI_targets[ ,15:21]- CI_targets[ ,1:7])/3.92
  SE_females <- (CI_targets[ ,22:28] - CI_targets[ ,8:14])/3.92
  SE <- cbind(SE_males, SE_females)
  colnames(SE) <- paste(colnames(Targets),"_SE", sep="")
  
  
#Define parameters to calibrate
Calibr_parameters <- Params[c("P.onset", "P.onset_low.risk","P.onset_age", "RR.onset_sex", "P.sympt.diag_LGBC", "P.sympt.diag_St1",
                              "P.sympt.diag_St2", "P.sympt.diag_St3", "P.sympt.diag_St4",
                              "P.sympt.diag_Age",
                               "shape.t.StI.StII", "shape.t.StII.StIII", "shape.t.StIII.StIV", 
                              "P.LGtoHGBC", "P.ungiag.dead",
                              "Mean.t.StI.StII", "Mean.t.StII.StIII", "Mean.t.StIII.StIV") ,]

v.param_names <- rownames(Calibr_parameters) # number of parameters to calibrate 
n_params <- length(v.param_names)

n_samples <- 50 # number of calibration runs


if(calibration_type == "Random"){
  # Goodness-of-fit outcomes
  m.GOF <-  matrix(nrow=n_samples, ncol = n_targets)
  colnames(m.GOF) <- paste0(c(v.target_names), "_fit")
  
  source("R calibration//R//Random//f.random.R")
  print(output)
} else if(calibration_type == "Bayes"){
  
  source("R calibration//R//MHA_functions.R")
  source("R calibration//R//MHA.R")
  
}
