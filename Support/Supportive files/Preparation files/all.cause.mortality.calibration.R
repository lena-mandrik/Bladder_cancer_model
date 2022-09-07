# Lena Mandrik

# 06.06.2022

setwd("C:\\Users\\cm1om\\Documents\\Research\\Cancer\\Bladder\\model data")

OC.m <-as.matrix(read.table("Data\\OC_mortality.txt",header=F))
OC.m <-OC.m[30:100,]
colnames(OC.m) <- c("OCM_males", "OCM_females")
rownames(OC.m) <- 30:100

population <-as.matrix(read.table("population2018.txt",header=T))

head(population)

# Assign RR for all-cause mortality for smokers and former smokers
m.RR.smoke <- matrix(data =c(1.2,1.14,1.26, 2.76, 2.71, 2.81), ncol=3, nrow=2) #replace with the value from the parameter sheet
colnames(m.RR.smoke) <- c("mean", "CI_low", "CI_up")
rownames(m.RR.smoke) <- c("RR_former_smoke", "RR_current_smoke")

#' Function to proceed the input data
#' @details
#' This function calculates individualised Bladder cancer and other cause mortality risk based on different factors 
#' @params
#' pop: population matrix containing individual level attributes
#' Param_sets: File containing sampled parameter sets
#' @return Vector of individualised cancer risks
RR_CurrSmoke <-1.2
RR_PastSmoke <- 2.76
  
pop <- population


#Calculate mean population attributes
smoke <- sum(pop[, "current_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])
psmoke <- sum(pop[, "past_smoke"] * pop[, "weighting"])/sum(pop[, "weighting"])

#Make a new temporary matrix where pop attributes match CRC risk categories
popA <- matrix(0, nrow = length(pop[,1]), ncol = 2)
colnames(popA) <- c("age_m", "age_f")

as.numeric(rownames(OC.m))

popA[,"age_m"] <- replace(popA[,"age_m"], pop[,"sex"] ==1, 1)
popA[,"age_f"] <- replace(popA[,"age_f"], pop[,"sex"] ==0, 1)


#Calculate individual OCM risk 
risk_OCM <- (RR_CurrSmoke ^ (pop[, "current_smoke"] - smoke)) *
  (RR_PastSmoke ^ (pop[, "past_smoke"] - psmoke))

pop[, "age_0"] <- 30
replace(risk_OCM, pop[, "sex"] ==1, OC.m[as.numeric(rownames(OC.m)) == pop[, "age_0"], ])



# Risk of OCM at age 30
RR_OCM <- m.RR.smoke["RR_former_smoke",1]^(population[,"past_smoke"] - OC.m[1,1])

# 01/06/2022

# The code is used to calibrate the all cause mortality among smokers, past smokers, and non smokers



setwd("C:\\Users\\cm1om\\Documents\\Research\\Cancer\\Bladder\\model data\\Data")

# Load all-cause mortality
All.cause.m <-as.matrix(read.table("All_cause_mortality_2021.txt",header=TRUE))

# Load BC mortality 
BC.m <-as.matrix(read.table("BC_mort.txt",header=TRUE))


# Assign RR of all-cause mortality for smokers and former smokers

RR_former_smoke <- 1.2
RR_current_smoke <- 2.76
RR_no_smoke_to_all_pop <- 1 #a parameter that will be calibrated


# in 2020 more men (15.3%) than women (13.7%) smoked in Great Britain (ONS)
smoke_men <- 0.153
smoke_women <- 0.137

# Create a table for the transformed mortality

Table.mortality <- matrix(nrow= nrow(All.cause.m), ncol = 13)

colnames(Table.mortality) <- c("all_cause_male", "all_cause_female", "P_smoke", "P_past_smoke", 
                               "P_non_smoke", "all_cause_male_no_smoke", "all_cause_male_smoke", "all_cause_male_past_smoke", 
                               "all_cause_female_no_smoke", "all_cause_female_smoke", "all_cause_female_past_smoke",
                               "all_cause_male_compare", "all_cause_female_compare")


head(Table.mortality)

# Replace with the average smoking rate among 30+ olds
Table.mortality <- Table.mortality[31:101,]

Table.mortality[,"P_smoke"] <- smoke_men


# replace the values from the all-cause mortality data set
Table.mortality[,1:5] <- All.cause.m

Table.mortality_new <- Table.mortality

RR_no_smoke_to_all_pop<-0.99

while(RR_no_smoke_to_all_pop > 0.2){

Table.mortality_new[,"all_cause_male_no_smoke"] <- Table.mortality_new[,"all_cause_male"]*RR_no_smoke_to_all_pop
Table.mortality_new[,"all_cause_male_smoke"] <- Table.mortality_new[,"all_cause_male_no_smoke"]*RR_current_smoke
Table.mortality_new[,"all_cause_male_past_smoke"] <- Table.mortality_new[,"all_cause_male_no_smoke"]*RR_former_smoke

Table.mortality_new[,"all_cause_female_no_smoke"] <- Table.mortality_new[,"all_cause_female"]*RR_no_smoke_to_all_pop
Table.mortality_new[,"all_cause_female_smoke"] <- Table.mortality_new[,"all_cause_female_no_smoke"]*RR_current_smoke
Table.mortality_new[,"all_cause_female_past_smoke"] <- Table.mortality_new[,"all_cause_female_no_smoke"]*RR_former_smoke

Table.mortality_new[,"all_cause_male_compare"] <- Table.mortality_new[,"all_cause_male_no_smoke"]*Table.mortality_new[,"P_non_smoke"]+Table.mortality_new[,"all_cause_male_smoke"]*Table.mortality_new[,"P_smoke"] + Table.mortality_new[,"all_cause_male_past_smoke"]*Table.mortality_new[,"P_past_smoke"] 

if(sum(Table.mortality_new[,"all_cause_male"]-Table.mortality_new[,"all_cause_male_compare"])^2 < 0.000001){print(RR_no_smoke_to_all_pop)}

RR_no_smoke_to_all_pop <- RR_no_smoke_to_all_pop - 0.0001
}


head(Table.mortality_new)
tail(Table.mortality_new)

# Check the best outcome
RR_no_smoke_to_all_pop <-0.652

Table.mortality_new[,"all_cause_male_no_smoke"] <- Table.mortality_new[,"all_cause_male"]*RR_no_smoke_to_all_pop
Table.mortality_new[,"all_cause_male_smoke"] <- Table.mortality_new[,"all_cause_male"]*RR_current_smoke
Table.mortality_new[,"all_cause_male_past_smoke"] <- Table.mortality_new[,"all_cause_male"]*RR_former_smoke

Table.mortality_new[,"all_cause_female_no_smoke"] <- Table.mortality_new[,"all_cause_female"]*RR_no_smoke_to_all_pop
Table.mortality_new[,"all_cause_female_smoke"] <- Table.mortality_new[,"all_cause_female"]*RR_current_smoke
Table.mortality_new[,"all_cause_female_past_smoke"] <- Table.mortality_new[,"all_cause_female"]*RR_former_smoke

Table.mortality_new[,"all_cause_male_compare"] <- Table.mortality_new[,"all_cause_male_no_smoke"]*Table.mortality_new[,"P_non_smoke"]+Table.mortality_new[,"all_cause_male_smoke"]*Table.mortality_new[,"P_smoke"] + Table.mortality_new[,"all_cause_male_past_smoke"]*Table.mortality_new[,"P_past_smoke"] 

sum(Table.mortality_new[,"all_cause_male"]-Table.mortality_new[,"all_cause_male_compare"])

Table.mortality_new[,c(1,12)]
