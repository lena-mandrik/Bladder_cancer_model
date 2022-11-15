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
rownames(OC_mort) <- c(paste("0",c(1:101), sep = ""), paste("1",c(1:101), sep = ""))

# Load bladder cancer mortality
BC.mort1 <- matrix(0, ncol = 1, nrow = 142)
BC.mort2 <- matrix(0, ncol = 60, nrow = 142)
BC.1.mort <- cbind(BC.mort1, as.matrix(S1), BC.mort2)
BC.2.mort <- cbind(BC.mort1, as.matrix(S2), BC.mort2)
BC.3.mort <- cbind(BC.mort1, as.matrix(S3), BC.mort2)
BC.4.mort <- cbind(BC.mort1, as.matrix(S4), BC.mort2)

colnames(BC.1.mort) <- colnames(BC.2.mort) <- colnames(BC.3.mort) <- colnames(BC.4.mort) <- c(0:70)
rownames(BC.1.mort) <- rownames(BC.2.mort) <- rownames(BC.3.mort) <- rownames(BC.4.mort) <- c(paste("0",c(30:100), sep = ""), paste("1",c(30:100), sep = ""))


#Generate parameter sets for PSA, adjusted for cycle length
Param_sets <- generate_parameters(Params, N_sets)
Param.names <- colnames(Param_sets)


#Load correlated parameter sets which are the calibrated parameters
# Corr_param_sets <- as.matrix(read.table("data/Corr_params.txt", header = TRUE))
#Add into parameter sets
# Param_sets[1:N_sets, 1:49] <- Corr_param_sets[1:N_sets, 1:49]

#Specify health states
states <- c("NoBC", "BC_LG", "BC_HG", "Death")
states_long <- c("NoBC", "BC_LG","St1_HG",  "DeathOC","St2_HG","St3_HG","St4_HG","DeathBC")
v.n <- c(1:4)

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
out_names <- c("TOTAL_COSTS", "BC_COSTS", "SCREEN_COSTS", "SURV_COSTS", "HARM_COSTS", 
               "QALYS", "LYS", "HG_St1_SYMPT", "HG_St2_SYMPT", "HG_St3_SYMPT", "HG_St4_SYMPT", "HG_St1_SCRN", "HG_St2_SCRN", 
               "HG_St3_SCRN", "HG_St4_SCRN", "LG_SCRN", "HG_St1_MORT", "HG_St2_MORT", "HG_St3_MORT", "HG_St4_MORT", "SCREEN_INVITE", "SCREEN_RESPOND", "SCREEN_FOLLOWED")


#create a matrix for all people in the dataset with the time they get each stage if they get invasive cancer
m.BC.T.to.Stage <- matrix(nrow = n.i, ncol = 3)
colnames(m.BC.T.to.Stage) <-c("T.onsetToStage2", "T.onsetToStage3", "T.onsetToStage4")
rownames(m.BC.T.to.Stage) <- 1:n.i


#Set up random number array for each individual
m.Rand <- generate_random()

