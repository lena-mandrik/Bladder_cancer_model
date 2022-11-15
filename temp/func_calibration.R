
#####################################################################

#####Function to extract N of people alive at each age group. 
### The function uses the output matrix TR as an input and sums number of either population alive by age
## TR (the input of the function) - is a matrix with transition rates (proportion of people in each health state which sums to 1)


# Calculate the size of the cohorts 
Cohort_m <- sum(population[,"sex"]==1)
Cohort_f <- sum(population[,"sex"]==0)

alive_m <- as.matrix(rowSums(TR_m[ , c(1:3,5:7)]))
alive_f <- as.matrix(rowSums(TR_f[ , c(1:3,5:7)]))

# Extracting from the lists 
list_results <- results_no_screen
l.Diag <- map(list_results, ~.x$m.Diag)
TR_m <- map(list_results, ~.x$TR_m)
TR_f <- map(list_results, ~.x$TR_f)

# Get pop alive pop and average from the lists
alive_m <- map(TR_m, function(x) as.matrix(Cohort_m*rowSums(x[ , c(1:3,5:7)])))
alive_f <- map(TR_f, function(x) as.matrix(Cohort_f*rowSums(x[ , c(1:3,5:7)])))
alive_m <- Reduce('+', alive_m) / (length(alive_m))
alive_f <- Reduce('+', alive_f) / (length(alive_f))


################################################################
# @ input: the output of run_simulation: m.Diag , sex 1 - males, 0 - females
# @ output: incidence counts 

sex=1

f.diag.outcomes <- function(m.Diag, sex){ #function for long format incidence 
  
  m.diagnosed <- subset(m.Diag, m.Diag[ , "BC_diag"] !=0 & m.Diag[ , "sex"] == sex, select = c("BC_diag", "age_diag", "stage_diag", "LG_BC_diag", "age_LG_BC_diag", "age_BC_death"), drop = FALSE)
  
  v.incidence <-cbind(c(30:100),matrix(0, nrow = 71, ncol = 5))
  
  for (i in 1: 71){
    v.incidence[i,6] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.incidence[i,1]),1])
    v.incidence[i,2] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.incidence[i,1] & m.diagnosed[,"stage_diag"]==1),1])
    v.incidence[i,3] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.CRCincidence[i,1] & m.diagnosed[,"stage_diag"]==2),1])
    v.incidence[i,4] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.CRCincidence[i,1] & m.diagnosed[,"stage_diag"]==3),1])
    v.incidence[i,5] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.CRCincidence[i,1] & m.diagnosed[,"stage_diag"]==4),1])
    
  } 
  v.incidence
}