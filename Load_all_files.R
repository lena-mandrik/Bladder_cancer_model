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

#Load probabilities of BC mortality data by age, sex, and stage

list_of_file_names <- list.files(path = "Data/survival/", recursive = TRUE,
                            pattern = ".txt$", 
                            full.names = TRUE) # Extract the names of the files


list_of_files <- lapply(list_of_file_names, read.table, sep = "\t", header =T) #Add the files to the list
names(list_of_files) <- tools::file_path_sans_ext(basename(list_of_file_names)) #Extract the files from the list


#Load other cause mortality data
OC_mort <- as.matrix(read.table("Data/OC_mortality.txt"))
OC_mort_m <- as.matrix(OC_mort[,1]); OC_mort_f <- as.matrix(OC_mort[,2]) #extract matrices for males and females
rownames(OC_mort_m) <- rownames(OC_mort_f) <- 1:101


