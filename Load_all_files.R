### Code to load up all functions and data files ###

#Read all functions from all scripts within the R folder

sapply(list.files("R",full.names=T), source)

#Load parameter data

Params <- as.matrix(read.table("Data/Parameters.txt", header = TRUE))

#Load population data 
# Population from HSE 2018 is the base case population as the last survey that contains the EQ5D values

population <- as.matrix(read.table("Data\\population2018.txt", header = TRUE))
population[,"PID"] <- c(1:length(population[,1]))
n.i <- nrow(population)  # number of simulated individuals
population[, "risk_age"] <- 100 # set population risk score age to 100 to ensure screening doesn't start early

#Adjust population age if cohort model is chosen
if(cohort ==1){
  population[, "age"] <- cohort_age
}

#Load probabilities of BC mortality data by age, sex, and stage

list_of_file_names <- list.files(path = "Data/survival/", recursive = TRUE,
                            pattern = ".txt$", 
                            full.names = TRUE)


list_of_files <- lapply(list_of_file_names, read.table, sep = "\t", header =T)

names(list_of_files) <- tools::file_path_sans_ext(basename(list_of_file_names))


