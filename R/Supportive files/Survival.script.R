## Lena Mandrik
### 18/07/2022

### Supportive function to convert the survival into the probability to die annually
### May be used to create the probabilistic survival within CI


#Load and process BC mortality data by age and sex for each stage
m.mortality <- as.matrix(read.table("Support Data/BC_survival.txt", header = TRUE)) # The age in the database reflects the start age up to thenext threshold

# Create empty matrices to gather survival data and to convert it to p. to survive till the next year
m.survival <- matrix(ncol= 10, nrow=100)
row.names(m.survival) <- 1:100
colnames(m.survival) <- 1:10
m.mortality_by.y <- m.survival

# Convert to the dataframe for the future operations
df.mortality <- as.data.frame(m.mortality)


# Function f.survival.calc use empty data frames as inputs and df.mortality reporting survival data for 1,5,10 years by age and sex
# Includes the type of the outcome to gather (stage and mean or CI, eg. S1 - stage 1) and sex (sex_0), where 0 - females and 1 - males
# It returns the matrix with the probability to die from BC for a specific stage and sex by year and age (up to the year 10)
# It is assumed that after Y10 the probability to die from BC is zero


f.survival.calc <- function(df.mortality, m.survival, m.mortality_by.y, outcome, sex_0){
  
  
  for(age_0 in seq(from = 25, to =75, by=10)){ # Age categories (the survival is assumed the same within the age categories)
    
    
    # Interpolate between 1 and 5 years
    
    df <- subset(df.mortality, sex==sex_0 & age < (age_0 + 10) & age >= age_0) #subset the dataframe based on sex and age
    
    survival_5 <- seq(df[1,outcome], df[2,outcome], length.out=5)
    survival_10 <- seq(df[2,outcome], df[3,outcome], length.out=6)
    
    v.survival5.10 <- c(survival_5, survival_10[2:6])
   
    # l.model <- lm(df[ ,outcome] ~ df[ ,2]) # Fit the linear model, year is the defining variable df[ ,2]
    # coeff.model = l.model$coefficients # extract the coefficients
    # b = coeff.model[[1]]
    # a = coeff.model[[2]]
    
    
    for(age_i in age_0:(age_0+9)){ #repeat the model for each age in the group
      
     # for(year in 1:10){
        
       # m.survival[age_i,year] <- b + a*year
      
      m.survival[age_i, ] <- v.survival5.10
        
      #  year = year +1 
      #}
      age_i = age_i +1
    }
 }
  
  # Calculate the probability to die during the first year based on the survival function, assuming the survival is 100% at year 0
  m.mortality_by.y[30:84, 1] <- (100 - m.survival[30:84, 1])/100
  
  for(i in 2:10){
    m.mortality_by.y[30:84, i] <- (m.survival[30:84, (i-1)] - m.survival[30:84, i])/m.survival[30:84, (i-1)]
  }
  
  # The survival data are up to the age 75. Assume that those above age 75 has the same probability to die as those aged 75 yo.
  m.mortality_by.y <- m.mortality_by.y[30:84,] 
  m.mortality_by.y <- as.data.frame(m.mortality_by.y)
  
  for(i in 1:16){
    m.mortality_by.y <- add_row(m.mortality_by.y, m.mortality_by.y[54,])
    
  }
  
  rownames(m.mortality_by.y) <- 30:100 #Name the rows to reflect the age
  
  m.mortality_by.y # Return the matrix with the probabilities
  
}

# Save the survival matrices for Stages by sex
S1_f <- f.survival.calc(df.mortality, m.survival, m.mortality_by.y, outcome = "S1", sex_0=0)
S2_f <- f.survival.calc(df.mortality, m.survival, m.mortality_by.y, outcome = "S2", sex_0=0)
S3_f <- f.survival.calc(df.mortality, m.survival, m.mortality_by.y, outcome = "S3", sex_0=0)
S4_f <- f.survival.calc(df.mortality, m.survival, m.mortality_by.y, outcome = "S4", sex_0=0)

write.table(S1_f,"Data/S1_f.txt", sep = "\t")
write.table(S2_f,"Data/S2_f.txt", sep = "\t")
write.table(S3_f,"Data/S3_f.txt", sep = "\t")
write.table(S4_f,"Data/S4_f.txt", sep = "\t")


S1_m <- f.survival.calc(df.mortality, m.survival, m.mortality_by.y, outcome = "S1", sex_0=1)
S2_m <- f.survival.calc(df.mortality, m.survival, m.mortality_by.y, outcome = "S2", sex_0=1)
S3_m <- f.survival.calc(df.mortality, m.survival, m.mortality_by.y, outcome = "S3", sex_0=1)
S4_m <- f.survival.calc(df.mortality, m.survival, m.mortality_by.y, outcome = "S4", sex_0=1)

write.table(S1_m,"Data/S1_m.txt", sep = "\t")
write.table(S2_m,"Data/S2_m.txt", sep = "\t")
write.table(S3_m,"Data/S3_m.txt", sep = "\t")
write.table(S4_m,"Data/S4_m.txt", sep = "\t")