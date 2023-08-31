# Plot the outcomes and targets by 5-year period

Calibr_parameters["P.onset",1]=5.5E-06 
Calibr_parameters["P.onset_low.risk",1]= 0.6376
Calibr_parameters["P.onset_age",1]=1.138
Calibr_parameters["RR.onset_sex",1]=3.64 

Calibr_parameters["P.sympt.diag_LGBC",1]=0.1
Calibr_parameters["P.sympt.diag_St1",1]=0.15
Calibr_parameters["P.sympt.diag_St2",1]=0.25
Calibr_parameters["P.sympt.diag_St3",1]=0.3 
Calibr_parameters["P.sympt.diag_St4",1]=0.34 
Calibr_parameters["P.sympt.diag_Age",1]= 0.91
Calibr_parameters["shape.t.StI.StII",1]=1
Calibr_parameters["shape.t.StII.StIII",1]= 1
Calibr_parameters["shape.t.StIII.StIV",1]= 1
Calibr_parameters["P.LGtoHGBC",1]= 2.555e-03
Calibr_parameters["P.ungiag.dead",1]=-0.1 #-0.0051 
Calibr_parameters["Mean.t.StI.StII",1]= 4.9
Calibr_parameters["Mean.t.StII.StIII",1]=4.2
Calibr_parameters["Mean.t.StIII.StIV",1]=4

source("R calibration\\Main_model_calibration.R") # run the model

# Load predictions of the epidemiologic data
output <- f.calibr.output(results_no_screen) # extract the outputs
Predict <- cbind(output$rate_outcomes_m, output$rate_outcomes_f)

# Load targets by 5-years
Targets_by5 <- read.table("R calibration/Targets/Calibration_targets_bladder_counts.txt", header = TRUE, row.names=1)

# Load pop by 5-years
England.pop_by5 <- as.matrix(read.table("Data/Pop.alive.by.age.txt", header = TRUE)) # Load targets by 5-years

# Scale to the population which is in 1,000
#Targets_by5[ ,1:7] <-Targets_by5[ ,1:7]*England.pop_by5[,"Males"]*1000
#Targets_by5[ ,8:14] <-Targets_by5[ ,8:14]*England.pop_by5[,"Females"]*1000


#Predict <- Predict[,2:12]

#Calculate the predictions by 5-year
outputs_by_5 <- matrix(0, ncol=(ncol(Predict)), nrow=nrow(Targets_by5))

# get the counter for matrix of data
n.out=1

for(i.row in 1:(nrow(Targets_by5)-1)){
  for(i.col in 1:ncol(Predict)){
    
    outputs_by_5[i.row,i.col] = mean(Predict[n.out:(n.out+5),i.col]) #
  }
  n.out=n.out+5
  
}

#name rows and columns
colnames(outputs_by_5) <- colnames(Predict)
rownames(outputs_by_5) <- rownames(Targets_by5)

outputs_by_5[ ,1:7] <- outputs_by_5[ ,1:7]*England.pop_by5[1:nrow(Targets_by5),"Males"]*1000
outputs_by_5[ ,8:14] <-outputs_by_5[ ,8:14]*England.pop_by5[1:nrow(Targets_by5),"Females"]*1000

# Calculate SE
CI_targets <- read.table("R calibration/Targets/CI_targets_bladder_counts.txt", header = TRUE) # Beta distribution for the CI # f.CI.calc(Targets)
#CI_targets <-CI_targets[1:14,]

#p_lower = Targets_by5 -1.96*sqrt(Targets_by5)
#colnames(p_lower) <- paste(colnames(Targets_by5),"_lb", sep="")

#p_higher = Targets_by5 +1.96*sqrt(Targets_by5)
#colnames(p_higher) <- paste(colnames(Targets_by5),"_ub", sep="")

#CI_targets <- cbind(p_lower, p_higher)

#CI_targets[CI_targets <0] <-0

# Format data for ggplot
CI_targets <- as.data.frame(CI_targets)
Age <- seq(from=30, to=90, by=5)
#Age <-30:90


Target_3.4.M =rowSums(Targets_by5[,c("inc3.male","inc4.male")])
Target_3.4.F =rowSums(Targets_by5[,c("inc3.fem","inc4.fem")])

Predict_3.4.M =rowSums(outputs_by_5[,c("BC_stage4_male","BC_stage3_male")])
Predict_3.4.F =rowSums(outputs_by_5[,c("BC_stage4_female","BC_stage3_female")])

Targets_names = colnames(Targets_by5)
Predict_names = colnames(outputs_by_5)

Targets <- cbind(Targets_by5, Target_3.4.M, Target_3.4.F); colnames(Targets) = c(Targets_names, "inc3.4.male", "inc3.4.female")
Predict <- cbind(outputs_by_5, Predict_3.4.M, Predict_3.4.F); colnames(Predict) =c(Predict_names, "BC_stage3.4_male", "BC_stage3.4_female")

Targets <- as.data.frame(cbind(Age, Targets))
Predict <- as.data.frame(cbind(Age, Predict))
CI_targets <-as.data.frame(cbind(Age, CI_targets))
# Plot the predictions vs targets
# Visualize 
# Total incidence rate for males
#Predict <- as.data.frame(cbind(Age, Predict[1:61,]))
#Targets <- as.data.frame(Targets[1:61,])

(target_i.males = ggplot(data = Targets,
                         aes(x = Age,
                             y = incidence.male)) +
    geom_errorbar(aes(ymin = CI_targets$incidence.male_lb, ymax = CI_targets$incidence.male_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_total_male), color = 'red') +
    
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = Bc_stage1_male), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_stage2_male), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_stage3_male), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_stage4_male), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 4, males"))

# Stage 3,4 combined

(target_i.3.4.males = ggplot(data = Targets,
                          aes(x = Age,
                             y = inc3.4.male)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_stage3.4_male), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_LG_male), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence LG males"))

# Total incidence rate for females

(target_i.females = ggplot(data = Targets,
                           aes(x = Age,
                               y = incidence.fem)) +
    geom_errorbar(aes(ymin = CI_targets$incidence.fem_lb, ymax = CI_targets$incidence.fem_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_total_female), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = Bc_stage1_female), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_stage2_female), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 2, females"))

# Stage 3,4 incidence females
(target_i.3.females = ggplot(data = Targets,
                               aes(x = Age,
                                   y = inc3.fem)) +
    geom_errorbar(aes(ymin = CI_targets$inc3.fem_lb, ymax = CI_targets$inc3.fem_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_stage3_female ), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 3, females"))


(target_i.3.4.females = ggplot(data = Targets,
                            aes(x = Age,
                               y = inc3.4.female)) +
   geom_point() +
   geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_stage3.4_female), color = 'red') +
   theme(
     panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 3,4, females"))

(target_i.4.females = ggplot(data = Targets,
                               aes(x = Age,
                                   y = inc4.fem)) +
    geom_errorbar(aes(ymin = CI_targets$inc4.fem_lb, ymax = CI_targets$inc4.fem_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_stage4_female ), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence Stage 4, females"))

# LG incidence females
(target_i.LG.female = ggplot(data = Targets,
                             aes(x = Age,
                                 y = LG.fem)) +
    geom_errorbar(aes(ymin = CI_targets$LG.fem_lb, ymax = CI_targets$LG.fem_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_LG_female), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "Incidence LG females"))


# mortality


# Mortality from BC males
(target_mort.male = ggplot(data = Targets,
                           aes(x = Age,
                               y = BC.death.male)) +
    geom_errorbar(aes(ymin = CI_targets$BC.death.male_lb, ymax = CI_targets$BC.death.male_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_death_male), color = 'red') +
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
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_death_female), color = 'red') +
    theme(
      panel.border = element_rect(fill = NA, color = 'black')
    ) + 
    labs(title = "Calibration target",
         x = "Age",
         y = "BC death, females"))

