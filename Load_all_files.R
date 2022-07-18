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

#Load and process BC mortality data by age and sex for each stage
m.mortality <- as.matrix(read.table("Data/BC_survival.txt", header = TRUE))

m.survival <- matrix(ncol= 10, nrow=100)
row.names(m.survival) <- 1:100
colnames(m.survival) <- 1:10

m.surv.m.S1 <- m.surv.m.S2 <-m.surv.m.S3 <- m.surv.m.S4 <- m.surv.f.S1 <- m.surv.f.S2 <-m.surv.f.S3 <- m.surv.f.S4 <- m.survival  

m.mortality_by.y <- m.survival

# logistic regression

df.mortality <- as.data.frame(m.mortality)


age_0 = 25
sex_0 = 0
#year = 0 

###########################
df <-df1


# df.mortality, m.survival, m.mortality_by.y, outcome, age_0, sex_0

f.survival.calc <- function(df.mortality, m.survival, m.mortality_by.y, outcome, age_0, sex_0){

for(age_0 in seq(from = 25, to =75, by=10)){
  
  
df <- subset(df.mortality, sex==sex_0 & age < (age_0 + 10) & age >= age_0)

l.model <- lm(df[ ,outcome] ~ df[ ,year])
coeff.model = l.model$coefficients
b = coeff.model[[1]]
a = coeff.model[[2]]


for(age_i in age_0:(age_0+9)){
  
for(year in 1:10){
    
   m.survival[age_i,year] <- b + a*year
    
    year = year +1 
  }
    age_i = age_i +1
}
}


m.mortality_by.y[30:84, 1] <- (100 - m.survival[30:84, 1])/100

for(i in 2:10){
  m.mortality_by.y[30:84, i] <- (m.survival[30:84, (i-1)] - m.survival[30:84, i])/m.survival[30:84, (i-1)]
}

m.mortality_by.y <- m.mortality_by.y[30:84,] 

m.mortality_by.y <- as.data.frame(m.mortality_by.y)

for(i in 1:16){
  m.mortality_by.y <- add_row(m.mortality_by.y, m.mortality_by.y[54,])
  
}
rownames(m.mortality_by.y) <- 30:100

m.mortality_by.y

}


f.survival.calc(df.mortality, m.survival, m.mortality_by.y, outcome = "S1", age_0=25, sex_0=0)
  