
m.sympt = matrix(ncol=4, nrow=n.t)

for(n in 30:99){
  m.sympt[n-29,1] = sum(m.Diag[,"HG_age_diag"]==n & m.M_8s[ ,n-28]==3)/sum(m.M_8s[ ,n-28]==3)
  m.sympt[n-29,2] = sum(m.Diag[,"HG_age_diag"]==n & m.M_8s[ ,n-28]==5)/sum(m.M_8s[ ,n-28]==5)
  m.sympt[n-29,3] = sum(m.Diag[,"HG_age_diag"]==n & m.M_8s[ ,n-28]==6)/sum(m.M_8s[ ,n-28]==6)
  m.sympt[n-29,4] = sum(m.Diag[,"HG_age_diag"]==n & m.M_8s[ ,n-28]==7)/sum(m.M_8s[ ,n-28]==7)
  
  
}

# For stage 1
# Set the sample size
sample_size <- 1000000
 
 # Generate the first 50% of values as a constant value of 3
constant_values <- rep(3, sample_size/2)

# Generate the next 10% of values as a uniform outside of the range 0-7, i.e assigned as 8-12
constant_values2 <- round(qbeta(runif(sample_size*0.1, 0.05, 0.95), 1, 1, lower.tail=F) * 8+5) 

# Generate the remaining 50% of values from a beta distribution
beta_values <- round(qbeta(runif(sample_size/2, 0.05, 0.95), 1, 1, lower.tail=F) *3+2)   # Scale and shift to range 1 to 7

 # Combine the constant values and beta values
population <- c(constant_values2, beta_values, constant_values)
     
 # Round all the values to the nearest whole number
population <- round(population)

# Ensure all values are non-negative
population <- pmax(population, 0)

# calculate mean
mean(population)
          