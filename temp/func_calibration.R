
#####################################################################

#####Function to extract N of people alive at each age group. 
### The function uses the output matrix TR as an input and sums number of either population alive by age
## TR (the input of the function) - is a matrix with transition rates (proportion of people in each health state which sums to 1)

# Function that extracts the outcputs from the list 
# The function return a list with incidence of HGBC, LGBC, and undiagnosed BC by age and proportion of LGBC progressed to HGBC

f.calibr.output <- function(list_results){
  
l.Diag <- map(list_results, ~.x$m.Diag)
TR_m <- map(list_results, ~.x$TR_m)
TR_f <- map(list_results, ~.x$TR_f)

# Get pop alive pop and average from the lists
alive_m <- map(TR_m, function(x) as.matrix(Cohort_m*rowSums(x[ , c(1:3,5:7)])))
alive_f <- map(TR_f, function(x) as.matrix(Cohort_f*rowSums(x[ , c(1:3,5:7)])))
alive_m <- Reduce('+', alive_m) / (length(alive_m))
alive_f <- Reduce('+', alive_f) / (length(alive_f))

# Calculate total incidence and by stages and average from the lists
outcomes_m <- map(l.Diag, f.diag.outcomes, 1)
outcomes_f <- map(l.Diag, f.diag.outcomes, 0)
outcomes_m <- (Reduce('+', outcomes_m) / (length(outcomes_m)))
outcomes_f <- (Reduce('+', outcomes_f) / (length(outcomes_f)))

# Calculate proportion of LGBC that progressed to HGBC
c.LG.to.HG <- map(l.Diag, f.LG.to.HG)
c.LG.to.HG <- (Reduce('+', c.LG.to.HG) / (length(c.LG.to.HG)))

# Calculate outcomes per population alive
rate_outcomes_m <- outcomes_m[,2:9]/as.vector(alive_m) # all outcomes, the column 1 is the age
rate_outcomes_f <- outcomes_f[,2:9]/as.vector(alive_f) # all outcomes, the column 1 is the age

rownames(rate_outcomes_m) <- rownames(rate_outcomes_f) <- 30:100

results <- list(rate_outcomes_m = rate_outcomes_m,
                rate_outcomes_f = rate_outcomes_f,
                c.LG.to.HG =c.LG.to.HG)

results
}


################################################################
# @ input: the output of run_simulation: m.Diag , sex 1 - males, 0 - females
# @ output:  counts 


# function to calculate the incidence of HG cancer, undiagnosed HG cancer, LG cancer, and mortality by age

f.diag.outcomes <- function(m.Diag, sex){  
  
# Subset those who were diagnosed for incidence and mortality
  m.diagnosed <- subset(m.Diag, (m.Diag[ , "BC_diag"] !=0 | m.Diag[ , "LG_BC_diag"] !=0)  & m.Diag[ , "sex"] == sex, select = c("BC_diag", "age_diag", "stage_diag", "LG_BC_diag", "age_LG_BC_diag", "age_BC_death"), drop = FALSE)
  v.outcomes <-cbind(c(30:100),matrix(0, nrow = 71, ncol = 7))
  
# Subset those who are undiagnosed for undiagnosed cancer prevalence
  m.undiagnosed <- subset(m.Diag, m.Diag[ , "yr_onset"] !=0 & m.Diag[ , "sex"] == sex, select = c("BC_diag", "age_diag", "stage_diag", "yr_diag", "yr_onset", "age_onset"), drop = FALSE)
  m.undiag <-cbind(c(30:100),matrix(0, nrow = 71, ncol = 1))
  
  
  for (i in 1: 71){
    # Calculate the incidence by age and sex for HG cancers
    v.outcomes[i,6] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.outcomes[i,1]),"BC_diag"])
    v.outcomes[i,2] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.outcomes[i,1] & m.diagnosed[,"stage_diag"]==1),"BC_diag"])
    v.outcomes[i,3] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.outcomes[i,1] & m.diagnosed[,"stage_diag"]==2),"BC_diag"])
    v.outcomes[i,4] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.outcomes[i,1] & m.diagnosed[,"stage_diag"]==3),"BC_diag"])
    v.outcomes[i,5] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==v.outcomes[i,1] & m.diagnosed[,"stage_diag"]==4),"BC_diag"])
    
    # Calculate the incidence of LG cancers by age and sex
    v.outcomes[i,7] <- sum(m.diagnosed[which(m.diagnosed[,"age_LG_BC_diag"]==v.outcomes[i,1]),"LG_BC_diag"])
    
    # Calculate the mortality from BC by age
    v.outcomes[i,8] <- sum(m.diagnosed[which(m.diagnosed[,"age_BC_death"]==v.outcomes[i,1]),"BC_diag"])
    
    # Calculate undiagnosed High grade cancer
    m.undiag[i,2] <- sum(m.undiagnosed[which(m.undiagnosed[,"age_onset"]<=m.undiag[i,1] & m.undiagnosed[,"age_diag"] >m.undiag[i,1]),"BC_diag"])
      } 
  
  m.outcomes <- cbind(v.outcomes, m.undiag[ ,2])
  colnames(m.outcomes) <- c("Age", "Bc_stage1", "BC_stage2", "BC_stage3", "BC_stage4", "BC_total", "BC_LG", "BC_death", "Undiag_HGBC")
                            
  m.outcomes
}

# Function to calculate the rate of LG BC progressed to HRBG from all LGBR

f.LG.to.HG <- function(m.Diag){  

  m.diagnosed <- subset(m.Diag, (m.Diag[ , "BC_diag"] !=0 | m.Diag[ , "LG_BC_diag"] !=0), select = c("BC_diag", "age_diag", "stage_diag", "LG_BC_diag", "age_LG_BC_diag", "age_BC_death"), drop = FALSE)
  
  # Calculate the proportion of LG cancers that progressed to HG cancers
  LG_to_HR_cancers <- sum(m.diagnosed[ ,"LG_BC_diag"]& m.diagnosed[ ,"BC_diag"])/sum(m.diagnosed[ ,"LG_BC_diag"])
  
  LG_to_HR_cancers
}


# function to transform Targets in the short format (by age group) to the long format (by each age)

f.targets.per.alive.long <- function(target.data, start_age, end_age, n_by) {
  
  output <- start_age:end_age
  Age = seq(start_age, end_age, by =n_by)
  
  target.data = cbind(Age, target.data)
  
  for(n_outcome in 2:ncol(target.data)){
    
    start <-0 
    
    for (i in 1:(nrow(target.data)-1)){
      
      extrapolation <- (approx(x = target.data[i:(i+1), "Age"], y = target.data[i:(i+1), n_outcome], n=n_by))[[2]]
      
      start <- c(start, extrapolation)
    }
    
    start <- c(start[-1],start[length(start)]) # there is an atrifact of extrapolation; need to figure out why
    
    output <- cbind(output,start)
    
  }
  colnames(output) <- colnames(target.data)
  
  output
  
}