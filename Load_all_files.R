### Code to load up all functions and data files ###
# Read all functions from all scripts within the R folder
sapply(list.files("R", full.names=T), source)

#################################################################################################################################
# Load parameter data
#################################################################################################################################
# For the kidney model, load the kidney parameters
if(disease=="kidney"| disease == "bladder_kidney"){Params_KC <- as.matrix(read.table("Data/kidney/Parameters.txt", header = TRUE, row.names=1))}

if(disease=="bladder"| disease == "bladder_kidney"){ Params_BC <- as.matrix(read.table("Data/bladder/Parameters.txt", header = TRUE, row.names=1))}

# Load general Parameters which are identical for bladder and kidney and will be saved in the Global space
Params <- as.matrix(read.table("Data/Parameters_general.txt", header = TRUE, row.names=1))

#################################################################################################################################
# Load population data
#################################################################################################################################
# Population from HSE 2018 is the base case population as the last survey that contains the EQ5D values
population <- as.matrix(read.table("Data/population2018.txt", header = TRUE))

n.i <- nrow(population)  # number of simulated individuals
# Load the states if the model is run as a multi-age cohort
#states_pop <- as.matrix(read.table("Data\\states_pop.txt", header = TRUE))
#states_pop <- states_pop[ ,1:9]
#colnames(states_pop) <- 1:9 # The state 9 is diagnosed with cancer

# if the population consists of current or former smokers, replace the population by sampling
if(char_pop =="c.smoke"){
  
  # Replace the population with the current smokers
  population <- population[population[ ,"current_smoke"]==1, ][sample(nrow(population[population[ ,"current_smoke"]==1, ]), n.i, replace = TRUE), ]

  } else if(char_pop =="all.smoke"){
  population <- population[population[ ,"current_smoke"]==1 | population[ ,"past_smoke"]==1, ][sample(nrow(population[population[ ,"current_smoke"]==1 | population[ ,"past_smoke"]==1, ]), n.i, replace = TRUE), ]
}
 
# if the population consists of current or men only, replace the population by sampling
if(sex =="men"){
  
  # Replace the population with the current smokers
  population <- population[population[ ,"sex"]==1, ][sample(nrow(population[population[ ,"sex"]==1, ]), n.i, replace = TRUE), ]
  
} else if(sex =="women"){
  population <- population[population[ ,"sex"]==0, ][sample(nrow(population[population[ ,"sex"]==0, ]), n.i, replace = TRUE), ]
}
 
population[,"PID"] <- c(1:length(population[,1]))


# Adjust population age if cohort model is chosen
if(cohort ==1){
  population[, "age"] <- cohort_age
}

#################################################################################################################################
# Load Survival
#################################################################################################################################

#Load probabilities of BC mortality data by age, sex, and stage (as separate files to avoid big matrices)
#(probability assumed to be 0 after 10 years, and 0 if undiagnosed apart from stage IV)
# Survival is saved into disease-specific environments

# For kidney and both diseases, load kidney survival data
if(disease=="kidney" | disease == "bladder_kidney"){
  
  list_of_file_names <- list.files(path = "Data/kidney/Survival/", recursive = TRUE,
                                                       pattern = ".txt$", 
                                                       full.names = TRUE) 
  
  # Process survival into the necessary format
  f.C.mort(list_of_file_names, set.envir=e.KC, "_KC")
  
} 

# For bladder and both diseases, load bladder survival data
if(disease=="bladder" | disease == "bladder_kidney"){
    
  list_of_file_names <- list.files(path = "Data/bladder/Survival/", recursive = TRUE,
                                   pattern = ".txt$", 
                                   full.names = TRUE) # Extract the names of the files
  
  # Process survival into the necessary format
  f.C.mort(list_of_file_names, set.envir=e.BC, "_BC")
}

#################################################################################################################################
#Load other cause mortality data
#################################################################################################################################

if(disease=="kidney"){
OC_mort <- as.matrix(read.table("Data/kidney/OC_mortality.txt"))
} else if(disease=="bladder"){
  OC_mort <- as.matrix(read.table("Data/bladder/OC_mortality.txt"))
}else if(disease=="bladder_kidney"){
  OC_mort <- as.matrix(read.table("Data/OC_mortality.txt"))
}
rownames(OC_mort) <- c(paste("0",c(1:101), sep = ""), paste("1",c(1:101), sep = ""))
##################################################################################################################################

#Load correlated parameter sets which are the calibrated parameters
# Corr_param_sets <- as.matrix(read.table("data/Corr_params.txt", header = TRUE))
#Add into parameter sets
# Param_sets[1:N_sets, 1:49] <- Corr_param_sets[1:N_sets, 1:49]

#Specify health states based on whether bladder, kidney, or both are modelled 


states <- c("NoBC", "BC_LG", "BC_HG", "KC", "Death")
states_long <- c("NoBC", "BC_LG","St1_HG",  "DeathOC","St2_HG","St3_HG","St4_HG","DeathBC")
v.n <- c(1:5)
v.n_long <- c(1:8)
n.s   <- length(states)  # the number of health states modelled as states
n.s_long   <- length(states_long)  # the number of all health states 


# State 1. no disease (ND)
# State 2. low grade (LG)
# state 3. high grade (HG) Stage 1 (HG_St1)
# State 4. Death OC
# State 5. high grade (HG) Stage 2 (HG_St2)
# State 6. high grade (HG) Stage 3 (HG_St3)
# State 7. high grade (HG) Stage 4 (HG_St4)
# State 8. Death from bladder cancer


#Set costs and QALYs discount weights
v.dwc <- 1 / (1 + d.c) ^ (0:n.t)   # calculate the cost discount weight based on the discount rate d.c    
v.dwe <- 1 / (1 + d.e) ^ (0:n.t)   # calculate the QALY discount weight based on the discount rate d.e



#Specify outcomes
out_names <- c("TOTAL_COSTS", "Cancer_COSTS", "DIAG_COSTS", "SCREEN_COSTS", 
               "QALYS", "LYS", "LG_SYMPT", "HG_St1_SYMPT", "HG_St2_SYMPT", "HG_St3_SYMPT", "HG_St4_SYMPT", "LG_SCRN", "HG_St1_SCRN", "HG_St2_SCRN", 
               "HG_St3_SCRN", "HG_St4_SCRN",  "HG_St1_MORT", "HG_St2_MORT", "HG_St3_MORT", "HG_St4_MORT", "Die_TURBT", 
               "Invite_DS","Respond_DS", "Positive_DS", "Respond_Cyst", "Diagnostic_Cyst", "TURBT",
               "FP", "FN")



#Generate parameter sets for PSA, adjusted for cycle length
# Generate parameters that are going to be saved in the Global space
Param_sets <- f.generate_parameters(Params, N_sets)
Param.names <- colnames(Param_sets)

# Generate disease-specific parameters
if(disease=="kidney"){
Param_sets_KC <- f.generate_parameters(Params_KC, N_sets)
Param.names_KC <- colnames(Param_sets_KC)
} else if(disease=="bladder"){
Param_sets_BC <- f.generate_parameters(Params_BC, N_sets)
Param.names_BC <- colnames(Param_sets_BC)
}else if(disease=="bladder_kidney"){
  Param_sets_BC <- f.generate_parameters(Params_BC, N_sets)
  Param_sets_KC <- f.generate_parameters(Params_KC, N_sets)
  Param.names_BC <- colnames(Param_sets_BC)
  Param.names_KC <- colnames(Param_sets_KC)}
