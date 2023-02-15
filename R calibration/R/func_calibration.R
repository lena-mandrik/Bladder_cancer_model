


#########################################################################################
# Function to format the calibration parameters using the fitted parameters from the random calibration as the starting set


# Load starting parameters and set up the step
#
f.load.start.param <- function(fitted_params, Calibr_parameters){
  
  # Load the calibrated parameters from the random calibration
  #fitted_params <- (read.table("R calibration/Outputs/parameters_fit.txt", header = F, row.names=1)) #exclude the parameter on RR death for no smokers, as it was calculated
  Calibr_parameters[1:nrow(fitted_params),1] <- fitted_params[,1]
  
  # Add the step to the matrix
  step_calibr <- matrix(0.3, ncol=1, nrow=nrow(Calibr_parameters))
  names.col_param <- colnames(Calibr_parameters)
  Calibr_parameters <- cbind(Calibr_parameters, step_calibr)
  colnames(Calibr_parameters) <- c(names.col_param, "step")
  
  # Update alpha and beta for parameters with beta distribution
  Calibr_parameters["P.onset", "Param1"]=((1-Calibr_parameters["P.onset", "Mean"])/(0.01^2)^2 -1/Calibr_parameters["P.onset", "Mean"])*Calibr_parameters["P.onset", "Mean"]^2
  Calibr_parameters["P.onset", "Param2"]=Calibr_parameters["P.onset", "Param1"]*(1/Calibr_parameters["P.onset", "Mean"]-1)
  
  Calibr_parameters["P.onset_low.risk", "Param1"]=((1-Calibr_parameters["P.onset_low.risk", "Mean"])/(0.01^2)^2 -1/Calibr_parameters["P.onset_low.risk", "Mean"])*Calibr_parameters["P.onset_low.risk", "Mean"]^2
  Calibr_parameters["P.onset_low.risk", "Param2"]=Calibr_parameters["P.onset_low.risk", "Param1"]*(1/Calibr_parameters["P.onset_low.risk", "Mean"]-1)
  
  Calibr_parameters["P.LGtoHGBC", "Param1"]=((1-Calibr_parameters["P.LGtoHGBC", "Mean"])/(0.01^2)^2 -1/Calibr_parameters["P.LGtoHGBC", "Mean"])*Calibr_parameters["P.LGtoHGBC", "Mean"]^2
  Calibr_parameters["P.LGtoHGBC", "Param2"]=Calibr_parameters["P.LGtoHGBC", "Param1"]*(1/Calibr_parameters["P.LGtoHGBC", "Mean"]-1)
  
  return(Calibr_parameters)
}




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
rate_outcomes_m <- outcomes_m[,2:ncol(outcomes_m)]/as.vector(alive_m) # all outcomes, the column 1 is the age
rate_outcomes_f <- outcomes_f[,2:ncol(outcomes_f)]/as.vector(alive_f) # all outcomes, the column 1 is the age

rownames(rate_outcomes_m) <- rownames(rate_outcomes_f) <- 30:100

colnames(rate_outcomes_m) <- paste0(colnames(rate_outcomes_m), "_male")
colnames(rate_outcomes_f) <- paste0(colnames(rate_outcomes_f), "_female")


results <- list(rate_outcomes_m = rate_outcomes_m,
                rate_outcomes_f = rate_outcomes_f,
                c.LG.to.HG =c.LG.to.HG)

results
}


################################################################
# @ input: the output of run_simulation: m.Diag , sex 1 - males, 0 - females
# @ output:  counts 


# function to calculate the incidence of HG cancer, LG cancer, and mortality by age

f.diag.outcomes <- function(m.Diag, sex){  
  
# Subset those who were diagnosed for incidence and mortality
  m.diagnosed <- subset(m.Diag, (m.Diag[ , "BC_diag"] !=0 | m.Diag[ , "LG_BC_diag"] !=0)  & m.Diag[ , "sex"] == sex, select = c("BC_diag", "age_diag", "stage_diag", "LG_BC_diag", "age_LG_BC_diag", "age_BC_death"), drop = FALSE)
  m.outcomes <-cbind(c(30:100),matrix(0, nrow = 71, ncol = 6))
  

  for (i in 1: 71){
    # Calculate the incidence by age and sex for HG cancers
    m.outcomes[i,5] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==m.outcomes[i,1]),"BC_diag"])
    m.outcomes[i,2] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==m.outcomes[i,1] & m.diagnosed[,"stage_diag"]==1),"BC_diag"])
    m.outcomes[i,3] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==m.outcomes[i,1] & m.diagnosed[,"stage_diag"]==2),"BC_diag"])
    m.outcomes[i,4] <- sum(m.diagnosed[which(m.diagnosed[,"age_diag"]==m.outcomes[i,1] & (m.diagnosed[,"stage_diag"]==3|m.diagnosed[,"stage_diag"]==4)),"BC_diag"])

    # Calculate the incidence of LG cancers by age and sex
    m.outcomes[i,6] <- sum(m.diagnosed[which(m.diagnosed[,"age_LG_BC_diag"]==m.outcomes[i,1]),"LG_BC_diag"])
    
    # Calculate the mortality from BC by age
    m.outcomes[i,7] <- sum(m.diagnosed[which(m.diagnosed[,"age_BC_death"]==m.outcomes[i,1]),"BC_diag"])
    
      } 
  
  colnames(m.outcomes) <- c("Age", "Bc_stage1", "BC_stage2", "BC_stage3.4", "BC_total", "BC_LG", "BC_death")
                            
  m.outcomes
}


################################################################
# @ input: the output of run_simulation: m.Diag , sex 1 - males, 0 - females
# @ output:  rates of undiagnosed HG cancer 

f.undiag.BC <- function(m.Diag, sex){  
  
 
  # Subset those who are undiagnosed for undiagnosed cancer prevalence
  m.undiagnosed <- subset(m.Diag, m.Diag[ , "yr_onset"] !=0 & m.Diag[ , "sex"] == sex, select = c("age_diag", "BC_diag", "age_onset"), drop = FALSE)
  m.undiag <-cbind(c(30:100),matrix(0, nrow = 71, ncol = 1))
  
  
  for (i in 1: 71){
       # Calculate undiagnosed High grade cancer
    m.undiag[i,2] <- sum(m.undiagnosed[which(m.undiagnosed[,"age_onset"]<=m.undiag[i,1] & m.undiagnosed[,"age_diag"] >m.undiag[i,1]),"BC_diag"])
  } 
  
  colnames(m.undiag) <- c("Age", "Undiag_HGBC")
  
  m.undiag
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



# Apply Wilson Score method to calculate CI

# Calculate the uncertainty 

f.CI.Wilson.calc <- function(Targets){
  
England.pop <- as.matrix(read.table("Data/Pop.alive.by.age.txt", header = TRUE)) # Load targets by 5-years
England.pop <-  f.targets.per.alive.long(England.pop, start_age=30, end_age =95, n_by=5) # the pop is in thousands

# for males
q =1 - as.matrix(Targets[,-c(1,7:11)])
dev.males = as.matrix(England.pop[,"Males"])*1000

denom = 2*(dev.males + 1.96^2)
nominator_l =(20 + 1.96^2 -1.96*sqrt(40*(q)+1.96^2))
p_lower_m =nominator_l /denom[col(nominator_l)]
nominator_h =(20 + 1.96^2 +1.96*sqrt(1.96^2 + 40*(q)))
p_higher_m = nominator_h /denom[col(nominator_h)]

# for females
q = 1 - (Targets[,c(7:11)])
dev.females = as.matrix(England.pop[,"Females"])*1000

denom = 2*(dev.females + 1.96^2)
nominator_l =(20 + 1.96^2 -1.96*sqrt(1.96^2 + 40*(q)))
p_lower_f =nominator_l /denom[col(nominator_l)]
nominator_h =(20 + 1.96^2 +1.96*sqrt(1.96^2 + 40*(q)))
p_higher_f = nominator_h /denom[col(nominator_h)]

colnames(p_higher_f) <- paste(colnames(p_higher_f),"_ub", sep="")
colnames(p_higher_m) <- paste(colnames(p_higher_m),"_ub", sep="")
colnames(p_lower_f) <- paste(colnames(p_lower_f),"_lb", sep="")
colnames(p_lower_m) <- paste(colnames(p_lower_m),"_lb", sep="")

output <- cbind(p_lower_m, p_higher_m, p_lower_f, p_higher_f)

output
}

# Function that calculates a standard CI based on the reported rates

f.CI.calc <- function(Targets){
  
  England.pop <- as.matrix(read.table("Data/Pop.alive.by.age.txt", header = TRUE)) # Load targets by 5-years
  England.pop <-  f.targets.per.alive.long(England.pop, start_age=30, end_age =95, n_by=5) # the pop is in thousands
  
  # for males
  p =as.matrix(Targets[,-c(1,8:13)])
  denom=sqrt(as.matrix(England.pop[,"Males"])*1000)
  
  nominator =1.96*sqrt(p*(1-p))
  bond = nominator /denom[col(nominator)]
  p_lower_m = p-bond
  p_higher_m = p+bond
  
  # for females
  p = as.matrix(Targets[,c(8:13)])
  denom=sqrt(as.matrix(England.pop[,"Females"])*1000)
  bond = nominator /denom[col(nominator)]
  nominator =1.96*sqrt(p*(1-p))
  p_lower_f = p-bond
  p_higher_f = p+bond
  
  colnames(p_higher_f) <- paste(colnames(p_higher_f),"_ub", sep="")
  colnames(p_higher_m) <- paste(colnames(p_higher_m),"_ub", sep="")
  colnames(p_lower_f) <- paste(colnames(p_lower_f),"_lb", sep="")
  colnames(p_lower_m) <- paste(colnames(p_lower_m),"_lb", sep="")
  
  output <- cbind(p_lower_m, p_higher_m, p_lower_f, p_higher_f)
  
  output[output <0] <-0
  
  output
}




# Function that visualizes the target by age
f.plot.target <- function(Targets, CI_targets){
  
  # Data sets
  Targets <- as.data.frame(Targets)
  CI_targets <- as.data.frame(CI_targets)
  
  # Total incidence rate for males
 
  (target_i.males = ggplot(data = Targets,
                      aes(x = Age,
                          y = incidence.rate.male)) +
    geom_errorbar(aes(ymin = CI_targets$incidence.rate.male_lb, ymax = CI_targets$incidence.rate.male_ub)) +
    geom_point() +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence, males"))
  
  # Stage 1 incidence males
  (target_i.1.males = ggplot(data = Targets,
                           aes(x = Age,
                               y = inc1.male)) +
      geom_errorbar(aes(ymin = CI_targets$inc1.male_lb, ymax = CI_targets$inc1.male_ub)) +
      geom_point() +
      theme(
        panel.border = element_rect(fill = NA, color = 'black')
      ) + 
      labs(title = "Calibration target",
           x = "Age",
           y = "Incidence Stage 1, males"))
  
  # Stage 2 incidence males
  (target_i.2.males = ggplot(data = Targets,
                             aes(x = Age,
                                 y = inc2.male)) +
      geom_errorbar(aes(ymin = CI_targets$inc2.male_lb, ymax = CI_targets$inc2.male_ub)) +
      geom_point() +
      theme(
        panel.border = element_rect(fill = NA, color = 'black')
      ) + 
      labs(title = "Calibration target",
           x = "Age",
           y = "Incidence Stage 2, males"))
  
  # Stage 3,4 incidence males
  (target_i.3.4.males = ggplot(data = Targets,
                             aes(x = Age,
                                 y = inc3.4.male)) +
      geom_errorbar(aes(ymin = CI_targets$inc3.4.male_lb, ymax = CI_targets$inc3.4.male_ub)) +
      geom_point() +
      theme(
        panel.border = element_rect(fill = NA, color = 'black')
      ) + 
      labs(title = "Calibration target",
           x = "Age",
           y = "Incidence Stage 3,4, males"))
  
  # LG incidence males
  (target_i.LG.male = ggplot(data = Targets,
                               aes(x = Age,
                                   y = LG.male)) +
      geom_errorbar(aes(ymin = CI_targets$LG.male_lb, ymax = CI_targets$LG.male_ub)) +
      geom_point() +
      theme(
        panel.border = element_rect(fill = NA, color = 'black')
      ) + 
      labs(title = "Calibration target",
           x = "Age",
           y = "Incidence LG males"))
  
  # Total incidence rate for females
  
  (target_i.females = ggplot(data = Targets,
                           aes(x = Age,
                               y = incidence.rate.fem)) +
      geom_errorbar(aes(ymin = CI_targets$incidence.rate.fem_lb, ymax = CI_targets$incidence.rate.fem_ub)) +
      geom_point() +
      theme(
        panel.border = element_rect(fill = NA, color = 'black')
      ) + 
      labs(title = "Calibration target",
           x = "Age",
           y = "Incidence, females"))
  
  # Stage 1 incidence females
  (target_i.1.females = ggplot(data = Targets,
                             aes(x = Age,
                                 y = inc1.fem)) +
      geom_errorbar(aes(ymin = CI_targets$inc1.fem_lb, ymax = CI_targets$inc1.fem_ub)) +
      geom_point() +
      theme(
        panel.border = element_rect(fill = NA, color = 'black')
      ) + 
      labs(title = "Calibration target",
           x = "Age",
           y = "Incidence Stage 1, females"))
  
  # Stage 2 incidence males
  (target_i.2.females = ggplot(data = Targets,
                             aes(x = Age,
                                 y = inc2.fem)) +
      geom_errorbar(aes(ymin = CI_targets$inc2.fem_lb, ymax = CI_targets$inc2.fem_ub)) +
      geom_point() +
      theme(
        panel.border = element_rect(fill = NA, color = 'black')
      ) + 
      labs(title = "Calibration target",
           x = "Age",
           y = "Incidence Stage 2, females"))
  
  # Stage 3,4 incidence females
  (target_i.3.4.females = ggplot(data = Targets,
                               aes(x = Age,
                                   y = inc3.4.fem)) +
      geom_errorbar(aes(ymin = CI_targets$inc3.4.fem_lb, ymax = CI_targets$inc3.4.fem_ub)) +
      geom_point() +
      theme(
        panel.border = element_rect(fill = NA, color = 'black')
      ) + 
      labs(title = "Calibration target",
           x = "Age",
           y = "Incidence Stage 3,4, females"))
  
  # LG incidence females
  (target_i.LG.female = ggplot(data = Targets,
                             aes(x = Age,
                                 y = LG.fem)) +
      geom_errorbar(aes(ymin = CI_targets$LG.fem_lb, ymax = CI_targets$LG.fem_ub)) +
      geom_point() +
      theme(
        panel.border = element_rect(fill = NA, color = 'black')
      ) + 
      labs(title = "Calibration target",
           x = "Age",
           y = "Incidence LG females"))
  
  print(list(target_i.males, target_i.1.males,target_i.2.males, target_i.3.4.males, target_i.LG.male,
        target_i.females, target_i.1.females, target_i.2.females, target_i.3.4.females, target_i.LG.female))
}

######################################################
# function to calculate GOF in the random calibration approach

f.GOF.calc <- function(run, m.GOF, Targets, SE, Predict){
  
  Targets <- Targets[ ,-1] #Get rid of the Age column, which the 1st one

  for(n in 1:ncol(Targets)){
    
    distribution <- dnorm(x = Targets[,n], mean = Predict[,n], sd = SE[,n], log = T)
    distribution[which(!is.finite(distribution))] <-0 # replace with zero infinite and undefined numbers
    m.GOF[run,n]<- sum(distribution)
    #m.GOF[run,n]<- sum(dnorm(x = Targets[,n], mean = Prediction[,n], sd = SE[,n], log = T))
  }
  m.GOF
}

