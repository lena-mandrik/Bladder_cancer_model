####
# Calculate the smoking cessation rate by 10y bands


# Set whether the population needs to be all population, current smokers, or both current and past smokers to be sampled before the model start
# Choose the population that is going to be simulated
char_pop <- "all" # c.smoke- current smokers, all.smoke- current and past smokers, all -the whole population
# Choose who from the simulated population is goign to receive the intervention
screen_elig <- "all" # c.smoke- current smokers, all.smoke- current and past smokers, all -the whole population
sex <- "all"

###Load up all the functions and all the data for use in the model
source("Load_all_files.R")

population <- as.matrix(read.table("Data/population2018.txt", header = TRUE))

age_band = seq(30,90,10)

output= matrix(nrow=length(age_band), ncol=4)
output[,1]=age_band
colnames(output) =c("age","past_smoke","current_smoke","no_smoke")

for(i in 1:length(age_band)){
  age = age_band[i]
  subset_pop <- population[population[,"age"]>=age & population[,"age"] <age+9, ]
  output[i,2] <- sum(subset_pop[,"past_smoke"])
  output[i,3] <- sum(subset_pop[,"current_smoke"])
  output[i,4] <- sum(subset_pop[,"no_smoke"])
  
}

output.cs = matrix(nrow=nrow(output), ncol=2)

output.cs = output[ ,2:4]/rowSums(output[ ,2:4])

output.cs=cbind(age_band, output.cs)
smoke_cs =smoke_ps =rep(0,6)

for(i in 1:(nrow(output.cs)-1)){
  smoke_cs[i]= ((output.cs[i,"current_smoke"]-output.cs[(i+1),"current_smoke"])/output.cs[i,"current_smoke"])/10
  smoke_ps[i]= ((output.cs[(i+1),"past_smoke"]-output.cs[i,"past_smoke"])/output.cs[(i+1),"past_smoke"])/10
  
}

write.csv(smoke_cs, file="Data\\smoke_cs.txt")
