## The script to analyse the output of the calibration and plot the predictions comparad to the targets
# Load GOF and parameters

#############################################################################################
  
fitted_param_distr <- (read.table("R calibration/Outputs/param.norm.txt", header = F, row.names=1))


for(i  in 1:nrow(Calibr_parameters)){
  p.mean =as.numeric(fitted_param_distr["mean" , i])
  p.sd = as.numeric(fitted_param_distr["sd" , i])
  Calibr_parameters[i,"Mean"] <-rnorm(1, p.mean, p.sd)
}



source("R calibration\\Main_model_calibration.R") # run the model


####################################

### For validation
# Calculate the size of the cohorts 
Cohort_m <- sum(population[,"sex"]==1)
Cohort_f <- sum(population[,"sex"]==0)

out.validation <- f.valid.output(results_no_screen)

time.to.symptoms <- out.validation$time.to.symptoms_diag
colnames(time.to.symptoms) <-c("Age", "time.to.symptoms")
time.to.symptoms <-as.data.frame(time.to.symptoms)

# Time to diagnosis
(valid.time.to.symptoms = ggplot(data = time.to.symptoms,
                         aes(x = Age,
                             y = time.to.symptoms)) + 
    geom_point() +
    geom_line(data= time.to.symptoms, aes(x = Age, y = time.to.symptoms), color = 'red') +
    labs(title = "Time to symptoms by age",
         x = "Age, years",
         y = "Time to symptoms among diagnosed, years"))

# Proportion of those who were diagnosed is low
p.diagnosed=out.validation$p.diagnosed

# With undiagnosed disease
d.undiag_m <- out.validation$undiag_m
#d.undiag_m <-d.undiag_m[,2]
Age= 30:100
d.undiag_m <- as.data.frame(cbind(Age, d.undiag_m))
colnames(d.undiag_m) <-c("Age", "undiag_m")

(plot.undiag_m = ggplot(data = d.undiag_m,
                                 aes(x = Age,
                                     y = undiag_m)) + 
    geom_point() +
    geom_line(data= d.undiag_m, aes(x = Age, y = undiag_m), color = 'red') +
    labs(title = "Undiagnosed HG BC by age, males",
         x = "Age, years",
         y = "Undiagnosed HG BC"))

d.undiag_f <- out.validation$undiag_f
#d.undiag_f <-d.undiag_f[,2]
Age= 30:100
d.undiag_f <- as.data.frame(cbind(Age, d.undiag_f))
colnames(d.undiag_f) <-c("Age", "undiag_f")

(plot.undiag_f = ggplot(data = d.undiag_f,
                        aes(x = Age,
                            y = undiag_f)) + 
    geom_point() +
    geom_line(data= d.undiag_f, aes(x = Age, y = undiag_f), color = 'red') +
    labs(title = "Undiagnosed HG BC by age, females",
         x = "Age, years",
         y = "Undiagnosed HG BC"))

# Prevalence of the disease
prev_m <- as.data.frame(out.validation$prev_m)*100000
HGBC <- rowSums(prev_m[ ,2:4])
Age= 30:100
prev_m <- as.data.frame(cbind(Age, prev_m[,1],HGBC))

colnames(prev_m) <-c("Age", "LG" ,"HG")

(plot.prev_m = ggplot(data = prev_m,
                        aes(x = Age,
                            y = LG)) + 
    geom_point() +
    geom_line(data= prev_m, aes(x = Age, y = HG), color = 'red') +
    labs(title = "Prevalence of LG and HG BC, males",
         x = "Age, years",
         y = "Prevalence BC"))



prev_f <- as.data.frame(out.validation$prev_f)*100000
HGBC <- rowSums(prev_f[ ,2:4])
Age= 30:100
prev_f <- as.data.frame(cbind(Age, prev_f[,1],HGBC))

colnames(prev_f) <-c("Age", "LG" ,"HG")

(plot.prev_f = ggplot(data = prev_f,
                      aes(x = Age,
                          y = LG)) + 
    geom_point() +
    geom_line(data= prev_f, aes(x = Age, y = HG), color = 'red') +
    labs(title = "Prevalence of LG and HG BC, females",
         x = "Age, years",
         y = "Prevalence BC"))

###################################################################
# Combine prevalence with incidence

## the code below is to have plots with outcomes as continious variables
output <- f.calibr.output(results_no_screen) # extract the outputs
Predict <- cbind(output$rate_outcomes_m, output$rate_outcomes_f)
names.predict <-colnames(Predict)

if(target_stat == "counts"){
  Predict <- sum_age_predictions(Predict, England.pop)
  colnames(Predict) <-names.predict
}

# log likelihood instead of sum of squared errors
#m.GOF <- f.GOF.calc(run, m.GOF, Targets, SE, Predict)
c.names <- colnames(Predict)
Predict <- cbind(as.numeric(rownames(Predict)),Predict)
colnames(Predict) <- colnames(Targets)



###################################
all.BC.m <- cbind(prev_m, Predict[ ,"BC_total_male"]*100000)
colnames(all.BC.m) <- c("Age", "LG_prev", "HG_prev", "HG_incidence")

(plot.BC_m = ggplot(data = all.BC.m,
                      aes(x = Age,
                          y = HG_incidence)) + 
    geom_point() +
    geom_line(data= all.BC.m, aes(x = Age, y = LG_prev), color = 'green') +
    geom_line(data= all.BC.m, aes(x = Age, y = HG_prev), color = 'red') +
    labs(title = "Prevalence and incidence of BC, males",
         x = "Age, years",
         y = "Prevalence and incidence BC"))


all.BC.f <- cbind(prev_f, Predict[ ,"BC_total_female"]*100000)
colnames(all.BC.f) <- c("Age", "LG_prev", "HG_prev", "HG_incidence")

(plot.BC_f = ggplot(data = all.BC.f,
                  aes(x = Age,
                      y = HG_incidence)) + 
    geom_point() +
    geom_line(data= all.BC.f, aes(x = Age, y = LG_prev), color = 'green') +
    geom_line(data= all.BC.f, aes(x = Age, y = HG_prev), color = 'red') +
    labs(title = "Prevalence and incidence of BC, females",
         x = "Age, years",
         y = "Prevalence and incidence BC"))
########################################################################

## If you want to plot the outcomes by 5-year slots, source the relevant plots file

if(target_stat == "rates"){
  
  # Extract the calibration targets and calculate the confidence interval
  # The incidence is already converted per 1 person alive
  Targets <- read.table("R calibration/Calibration_targets.txt", header = TRUE, row.names=1) # Load targets by 5-years
  Targets <- f.targets.per.alive.long(Targets, start_age=30, end_age =100, n_by=5) #Linearly extrapolate between the 5-year bands
  
  n_targets <- ncol(Targets)  # Exclude Age column
  v.target_names <- colnames(Targets) # Exclude Age column
  #CI_targets <- f.CI.calc(Targets)
  
  # Calculate SE from the CI
  SE_males <- (CI_targets[ ,15:21] - CI_targets[ ,1:7])/3.92
  SE_females <- (CI_targets[ ,22:28] - CI_targets[ ,8:14])/3.92
  SE <- cbind(SE_males, SE_females)
  colnames(SE) <- paste(colnames(Targets),"_SE", sep="")
  
  Age =30:84
  
  # Exclude pop age >95
  
  Targets <- Targets[1:55,]
  CI_targets <- CI_targets[1:55,]
  Predict <-Predict[1:55,]
  colnames(Predict) <- colnames(Targets[1:15])
} else if(target_stat == "counts"){
  
  Age =seq(30,90, by=5)
}

# Data sets
Targets <- as.data.frame(Targets)
CI_targets <- as.data.frame(CI_targets)



Targets <- cbind(Age, Targets)
CI_targets <-cbind(Age, CI_targets)



#Visualize by 5-year periods

#source("R calibration/R/plot_by_5y.R")
# Add the stage 3,4 combined
Targets_st3.4_m <- Targets[, "inc3.male"]+Targets[, "inc4.male"]
Targets_st3.4_f <- Targets[, "inc3.fem"]+Targets[, "inc4.fem"]
CI_lb_st3.4_m <- CI_targets[ ,"inc3.male_lb"]+CI_targets[ ,"inc4.male_lb"]
CI_ub_st3.4_m <- CI_targets[ ,"inc3.male_ub"]+CI_targets[ ,"inc4.male_ub"]
CI_lb_st3.4_f <- CI_targets[ ,"inc3.fem_lb"]+CI_targets[ ,"inc4.fem_lb"]
CI_ub_st3.4_f <- CI_targets[ ,"inc3.fem_ub"]+CI_targets[ ,"inc4.fem_ub"]
Predict_st3.4_m <- Predict[ ,"inc3.male"]+Predict[ ,"inc4.male"]
Predict_st3.4_f <- Predict[ ,"inc3.fem"]+Predict[ ,"inc4.fem"]

Targets <- cbind(Targets, Targets_st3.4_m, Targets_st3.4_f)
CI_targets <-cbind(CI_targets, CI_lb_st3.4_m, CI_ub_st3.4_m, CI_lb_st3.4_f, CI_ub_st3.4_f)
Predict <- cbind(Predict, Predict_st3.4_m, Predict_st3.4_f)

# visualise by year

# LG incidence females
(target_i.LG.female = ggplot(data = Targets,
                             aes(x = Age,
                                 y = LG.fem)) +
    geom_errorbar(aes(ymin = CI_targets$LG.fem_lb, ymax = CI_targets$LG.fem_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = LG.fem), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence LG females"))

# Total incidence rate for females

(target_i.females = ggplot(data = Targets,
                           aes(x = Age,
                               y = incidence.rate.fem )) +
    geom_errorbar(aes(ymin = CI_targets$incidence.rate.fem_lb, ymax = CI_targets$incidence.rate.fem_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = incidence.rate.fem), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = inc1.fem), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = inc2.fem), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 2, females"))

# Stage 3 incidence females
(target_i.3.females = ggplot(data = Targets,
                               aes(x = Age,
                                   y = inc3.fem)) +
    geom_errorbar(aes(ymin = CI_targets$inc3.fem_lb, ymax = CI_targets$inc3.fem_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = inc3.fem ), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 3, females"))

# Stage 4 incidence females
(target_i.4.females = ggplot(data = Targets,
                             aes(x = Age,
                                 y = inc4.fem)) +
    geom_errorbar(aes(ymin = CI_targets$inc4.fem_lb, ymax = CI_targets$inc4.fem_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = inc4.fem ), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 4, females"))


# Stage 3,4 combined incidence females
(target_i.3.4.females = ggplot(data = Targets,
                             aes(x = Age,
                                 y = Targets_st3.4_f)) +
    geom_errorbar(aes(ymin = CI_targets$CI_lb_st3.4_f, ymax = CI_targets$CI_ub_st3.4_f)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = Predict_st3.4_f), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 3,4, females"))


# Stage 3,4 combined incidence males
(target_i.3.4.males = ggplot(data = Targets,
                               aes(x = Age,
                                   y = Targets_st3.4_m)) +
    geom_errorbar(aes(ymin = CI_targets$CI_lb_st3.4_m, ymax = CI_targets$CI_ub_st3.4_m)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = Predict_st3.4_m), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 3,4, males"))


# Total incidence rate for males

# LG incidence males
(target_i.LG.male = ggplot(data = Targets,
                           aes(x = Age,
                               y = LG.male)) +
    geom_errorbar(aes(ymin = CI_targets$LG.male_lb, ymax = CI_targets$LG.male_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = LG.male), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence LG males"))

(target_i.males = ggplot(data = Targets,
                         aes(x = Age,
                             y = incidence.rate.male)) +
    geom_errorbar(aes(ymin = CI_targets$incidence.rate.male_lb, ymax = CI_targets$incidence.rate.male_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = incidence.rate.male), color = 'red') +
    
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
    geom_point()  +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = inc1.male), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = inc2.male), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 2, males"))

# Stage 3 incidence males
(target_i.3.males = ggplot(data = Targets,
                           aes(x = Age,
                               y = inc3.male)) +
    geom_errorbar(aes(ymin = CI_targets$inc3.male_lb, ymax = CI_targets$inc3.male_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = inc3.male), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 3, males"))

# Stage 4 incidence males
(target_i.4.males = ggplot(data = Targets,
                           aes(x = Age,
                               y = inc4.male)) +
    geom_errorbar(aes(ymin = CI_targets$inc4.male_lb, ymax = CI_targets$inc4.male_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = inc4.male), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 4, males"))



# mortality


# Mortality from BC males
(target_mort.male = ggplot(data = Targets,
                             aes(x = Age,
                                 y = BC.death.male)) +
    geom_errorbar(aes(ymin = CI_targets$BC.death.male_lb, ymax = CI_targets$BC.death.male_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC.death.male), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "BC death, males"))

# Mortality from BC females
(target_mort.female = ggplot(data = Targets,
                           aes(x = Age,
                               y = BC.death.female)) +
    geom_errorbar(aes(ymin = CI_targets$BC.death.female_lb, ymax = CI_targets$BC.death.female_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC.death.female), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "BC death, females"))


