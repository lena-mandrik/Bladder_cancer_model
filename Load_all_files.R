### Code to load up all functions and data files ###

# Read all functions from all scripts within the R folder

sapply(list.files("R",full.names=T), source)

# Load parameter data
Params <- as.matrix(read.table("Data/Parameters.txt", header = TRUE, row.names=1))

# Load population data 
# Population from HSE 2018 is the base case population as the last survey that contains the EQ5D values

population <- as.matrix(read.table("Data\\population2018.txt", header = TRUE))
population[,"PID"] <- c(1:length(population[,1]))
n.i <- nrow(population)  # number of simulated individuals
population[, "risk_age"] <- 100 # set population risk score age to 100 to ensure screening doesn't start early

# Adjust population age if cohort model is chosen
if(cohort ==1){
  population[, "age"] <- cohort_age
}

#Load probabilities of BC mortality data by age, sex, and stage (as separate files to avoid big matrices)
#(probability assumed to be 0 after 10 years, and 0 if undiagnosed apart from stage IV)

list_of_file_names <- list.files(path = "Data/survival/", recursive = TRUE,
                            pattern = ".txt$", 
                            full.names = TRUE) # Extract the names of the files


list_of_files <- lapply(list_of_file_names, read.table, sep = "\t", header =T) #Add the files to the list
names(list_of_files) <- tools::file_path_sans_ext(basename(list_of_file_names))  # Save the names
list2env(list_of_files,envir=.GlobalEnv) #Extract the files from the list

#Load other cause mortality data
OC_mort <- as.matrix(read.table("Data/OC_mortality.txt"))
OC_mort_m <- as.matrix(OC_mort[,1]); OC_mort_f <- as.matrix(OC_mort[,2]) #extract matrices for males and females
rownames(OC_mort_m) <- rownames(OC_mort_f) <- 1:101

#Generate parameter sets for PSA, adjusted for cycle length
Param_sets <- generate_parameters(Params, N_sets)
Param.names <- colnames(Param_sets)


#Load correlated parameter sets which are the calibrated parameters
# Corr_param_sets <- as.matrix(read.table("data/Corr_params.txt", header = TRUE))
#Add into parameter sets
# Param_sets[1:N_sets, 1:49] <- Corr_param_sets[1:N_sets, 1:49]

#Specify health states
states <- c("ND", "LG", "HG_St1", "HG_St2", "HG_St3", "HG_St4", "Diagnosed", "DeathBC", "DeathOC")
v.n <- c(1:9)
n.s   <- length(states)  # the number of health states
# State 1. no disease (ND)
# State 2. low grade (LG)
# state 3. high grade (HG) Stage 1 (HG_St1)
# State 4. high grade (HG) Stage 2 (HG_St2)
# State 5. high grade (HG) Stage 3 (HG_St3)
# State 6. high grade (HG) Stage 4 (HG_St4)
# State 7. Diagnosed

# State 8. DeathBC (Death from bladder cancer)
# State 9. DeathOC (Death from causes other than bladder cancer)

#Set costs and QALYs discount weights
v.dwc <- 1 / (1 + d.c) ^ (0:n.t)   # calculate the cost discount weight based on the discount rate d.c    
v.dwe <- 1 / (1 + d.e) ^ (0:n.t)   # calculate the QALY discount weight based on the discount rate d.e

#Specify outcomes
out_names <- c("TOTAL_COSTS", "BC_COSTS", "SCREEN_COSTS", "SURV_COSTS", "HARM_COSTS", 
               "QALYS", "LYS", "HG_St1_SYMPT", "HG_St2_SYMPT", "HG_St3_SYMPT", "HG_St4_SYMPT", "HG_St1_SCRN", "HG_St2_SCRN", 
               "HG_St3_SCRN", "HG_St4_SCRN", "LG_SCRN", "HG_St1_MORT", "HG_St2_MORT", "HG_St3_MORT", "HG_St4_MORT", "SCREEN_INVITE", "SCREEN_RESPOND", "SCREEN_FOLLOWED")
