# Lena Mandrik 12.07.2022
# This is the main script
### Run this script to run whole model
### Remove everything in workspace 
rm(list = ls(all = TRUE))

######################
###############

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
run_mode <- "Deterministic" # Available modes include "Testing" (returns all matrices), "Deterministic" (m.Out only), "PSA" (m.Out only)
cohort <- 1 # 1 = all individuals start model at same age (cohort), 0 = individuals start in model at true (HSE) age
disease <- "kidney" #"bladder" #"kidney" "bladder_kidney"

if(cohort ==1){
  cohort_age <- 30 #select starting age of cohort (hash out or set to anything if not using cohort)
  n.t   <- 100-cohort_age # The number of cycles to run 
} else if(cohort ==0){
  #set the min and max ages of the cohort if multi-age pop is simulated
  min_age<- 60
  max_age<- 80
  cohort_age<- 100

}


n.loops <- 2 # The number of model loops/PSA loops to run 
cl <- 1  # The cycle length (years) 
n.t   <- if(cohort==1){100-cohort_age}else{100-min_age}  # The number of cycles to run 
d.c <- 0.035 # The discount rate for costs
d.e <- 0.035 # The discount rate for effects
N_sets <- if(run_mode =="PSA"){n.loops}else{1} #Number of parameter sets required (for PSA). Minimum value = 1 (mean parameter values); maximum value = 1651 (number of calibrated correlated param sets)


# Set whether the population needs to be all population, current smokers, or both current and past smokers to be sampled before the model start
char_pop <- "all.smoke" # c.smoke- current smokers, all.smoke- current and past smokers, all -the whole population
screen_elig <- "all.smoke" # c.smoke- current smokers, all.smoke- current and past smokers, all -the whole population
sex <- "all"

###Load up all the functions and all the data for use in the model
source("Load_all_files.R")
nsample <- if(cohort==1){n.i} else{5000} #define a sample size for population to run in each loop 

###Load up all the functions and all the data for use in the model
###Set up results collection
results_no_screen <- results_screen_70 <- results_screen_75 <-results_screen_68 <- 
  results_screen_66 <-results_screen_64 <- results_screen_62 <-results_screen_60 <-list()

# create a list to collect PSA outcomes if the model runs in the PSA mode
if(run_mode == "PSA"){
  n_out = 1 # Add the number of sub populations for which the model is predicting the outcomes
  PSA_results_no_screen <- matrix(0, nrow = n_out * length(out_names), ncol = n.loops)
  rownames(PSA_results_no_screen) <- rep(out_names, n_out)
  PSA_results_screen_70 <- PSA_results_screen_72 <- PSA_results_screen_30 <- 
    PSA_results_screen_60_5t.1 <- PSA_results_screen_60_5t.2 <- PSA_results_screen_65_5t.1 <- 
    PSA_results_screen_65_5t.2 <- PSA_results_no_screen
}


#weight_outcome <- if(cohort==1){population[, "weighting"]} else{rep(1, nsample)} #Assign weights only if the cohort is modelled 


### run the model:
system.time(for(iter in 1:n.loops) {
  
#Select appropriate parameter set to use
p.set <- ifelse(run_mode =="PSA", iter, 1) 
#p.set=1

#Set up the model parameters according to current parameter set
f.set_parameters(p.set)

# Allocate the time to stage at diagnosis for each person in HSE
m.BC.T.to.Stage <- f.stage(Mean.t.StI.StII, shape.t.StI.StII, Mean.t.StII.StIII, shape.t.StII.StIII, 
                           Mean.t.StIII.StIV, shape.t.StIII.StIV, nsample)

#calculates relative and absolute risks of cancer onset
#pop <- f.risk.calc(population) 
pop <- set_population(population) 


pop_ns <- pop_sc_75 <- pop_sc_70 <- 
  pop_sc_68 <- pop_sc_66 <-pop_sc_64 <- pop_sc_62 <- pop_sc_60 <-pop #copy the population to ensure that individual characteristics remain as in the start

#Set up random number array for each individual
m.Rand <- f.generate_random()

# Update the smoking status for the previous smokers with the current smokers
pop <- f.smoke.initial(pop, m.Rand)

# temp code to check
#test_accuracy[c(2,3,5:7),1] <- 1
#diag1_accuracy[c(2,3,5:7),1] <- 1
#diag2_accuracy[c(2,3,5:7),1] <- 1
#test_accuracy[1,1] <- 0
#diag1_accuracy[1,1] <- 0
#diag2_accuracy[1,1] <- 0
#Mort.TURBT <-0

P.onset =5.5E-06 
P.onset_low.risk =0.638
P.onset_age =1.138
RR.onset_sex =3.64
P.sympt.diag_LGBC =0.102
P.sympt.diag_St1 =0.15
P.sympt.diag_St2=0.25
P.sympt.diag_St3=0.33
P.sympt.diag_St4=0.34
P.sympt.diag_Age =0.91
shape.t.StI.StII = 1
shape.t.StII.StIII =1
shape.t.StIII.StIV=1
P.LGtoHGBC=0.00255
P.ungiag.dead =-0.0051 
Mean.t.StI.StII=4.9
Mean.t.StII.StIII=4.2
Mean.t.StIII.StIV=4



# Run the model for no screening pop
DS_screen =0 #Set whether the screening with dipstick happens, 0 - no, 1- yes
DS_age = 100 #Set the age of the dipstick if screening happens, set to any if no screening
results_no_screen[[iter]] = Simulate_NHD(nsample, n.t, pop_ns, m.BC.T.to.Stage)

# Run the model for screening pop
DS_screen =1 #Set whether the screening with dipstick happens, 0 - no, 1- yes
DS_age = 75 #Set the age of the dipstick if screening happens, set to zero if no screening
DS_round = 1 #' DS_round - number of the screening rounds
DS_freq =0   #' DS_freq - frequency of the screening rounds (either 1 (annual) or 2 (biennial))


results_screen_75[[iter]] = Simulate_NHD(nsample, n.t, pop_sc_75, m.BC.T.to.Stage)

#DS_screen =1 #Set whether the screening with dipstick happens, 0 - no, 1- yes
#DS_age = 72 #Set the age of the dipstick if screening happens, set to zero if no screening
#DS_round = 1 #' DS_round - number of the screening rounds
#DS_freq =0   #' DS_freq - frequency of the screening rounds (either 1 (annual) or 2 (biennial))

#results_screen_72[[iter]] = Simulate_NHD(nsample, n.t, pop_sc_72, m.BC.T.to.Stage)

DS_screen =1 #Set whether the screening with dipstick happens, 0 - no, 1- yes
DS_age = 70 #Set the age of the dipstick if screening happens, set to zero if no screening
DS_round = 1 #' DS_round - number of the screening rounds
DS_freq =0   #' DS_freq - frequency of the screening rounds (either 1 (annual) or 2 (biennial))

results_screen_70[[iter]] = Simulate_NHD(nsample, n.t, pop_sc_70, m.BC.T.to.Stage)

DS_screen =1 #Set whether the screening with dipstick happens, 0 - no, 1- yes
DS_age = 68 #Set the age of the dipstick if screening happens, set to zero if no screening
DS_round = 1 #' DS_round - number of the screening rounds
DS_freq =0   #' DS_freq - frequency of the screening rounds (either 1 (annual) or 2 (biennial))

results_screen_68[[iter]] = Simulate_NHD(nsample, n.t, pop_sc_68, m.BC.T.to.Stage)

DS_screen =1 #Set whether the screening with dipstick happens, 0 - no, 1- yes
DS_age = 66 #Set the age of the dipstick if screening happens, set to zero if no screening
DS_round = 1 #' DS_round - number of the screening rounds
DS_freq =0   #' DS_freq - frequency of the screening rounds (either 1 (annual) or 2 (biennial))

results_screen_66[[iter]] = Simulate_NHD(nsample, n.t, pop_sc_66, m.BC.T.to.Stage)

DS_screen =1 #Set whether the screening with dipstick happens, 0 - no, 1- yes
DS_age = 64 #Set the age of the dipstick if screening happens, set to zero if no screening
DS_round = 1 #' DS_round - number of the screening rounds
DS_freq =0   #' DS_freq - frequency of the screening rounds (either 1 (annual) or 2 (biennial))

results_screen_64[[iter]] = Simulate_NHD(nsample, n.t, pop_sc_64, m.BC.T.to.Stage)

#DS_screen =1 #Set whether the screening with dipstick happens, 0 - no, 1- yes
#DS_age = 62 #Set the age of the dipstick if screening happens, set to zero if no screening
#DS_round = 1 #' DS_round - number of the screening rounds
#DS_freq =0   #' DS_freq - frequency of the screening rounds (either 1 (annual) or 2 (biennial))

#results_screen_62[[iter]] = Simulate_NHD(nsample, n.t, pop_sc_62, m.BC.T.to.Stage)

#DS_screen =1 #Set whether the screening with dipstick happens, 0 - no, 1- yes
#DS_age = 60 #Set the age of the dipstick if screening happens, set to zero if no screening
#DS_round = 1 #' DS_round - number of the screening rounds
#DS_freq =0   #' DS_freq - frequency of the screening rounds (either 1 (annual) or 2 (biennial))

#results_screen_60[[iter]] = Simulate_NHD(nsample, n.t, pop_sc_60, m.BC.T.to.Stage)

cat('\r', paste(round(iter/n.loops*100), "% done", sep = " "))


if(run_mode == "PSA"){
  PSA_results_no_screen[, iter] <- rowSums(results_no_screen[[iter]])
  PSA_results_screen_70[, iter] <- rowSums(results_screen_70[[iter]])
  PSA_results_screen_75[, iter] <- rowSums(results_screen_75[[iter]])
  PSA_results_screen_68[, iter] <- rowSums(results_screen_68[[iter]])

}

})



############################################
#' Analyse the results depending whether the runds were deterministic or probabilistic 
#' Analyse the results
if(run_mode == "Deterministic"){
  
  results_no_screen_cum <- (process_DA_results(results_no_screen, "results_NoScreen.txt"))
  results_screen_75_cum <- (process_DA_results(results_screen_75, "results_screen_70.txt"))
  results_screen_70_cum <- (process_DA_results(results_screen_70, "results_screen_72.txt"))
  results_screen_68_cum <- process_DA_results(results_screen_68, "results_screen_68.txt")
 results_screen_66_cum <- process_DA_results(results_screen_66, "results_screen_66.txt")
  results_screen_64_cum <- process_DA_results(results_screen_64, "results_screen_64.txt")
  #results_screen_62_cum <- (process_DA_results(results_screen_62, "results_screen_62.txt"))
 # results_screen_60_cum <- (process_DA_results(results_screen_60, "results_screen_60.txt"))
  
}

if(run_mode == "PSA"){
  
  results_no_screen_PSA_cum  <- process_PSA_results(results_no_screen, "results_no_screen_PSA.txt")
  results_no_screen_cum <- rowMeans(PSA_results_no_screen)
  results_screen_70_cum  <- rowMeans(PSA_results_screen_70)
  results_screen_75_cum  <- rowMeans(PSA_results_screen_75)
  PSA_results <- cbind(results_no_screen_mean, results_screen_70_mean, results_screen_75_mean)
}


########################
