## The script to analyse the output of the calibration and plot the predictions comparad to the targets
# Load GOF and parameters

#outputs <- read.table("R calibration/Outputs/m.GOF_21.12.22.txt", header = TRUE, row.names=1)

# Check fitted parameters from the previous calibration
fitted_params <- (read.table("R calibration/Outputs/parameters_fit.txt", header = F, row.names=1)) 
Calibr_parameters[ ,1] <- fitted_params[,1]

#Calibrated_parameters <- t(outputs[outputs[ ,"Sum"] ==max(outputs[ ,"Sum"]),14:25])

#Calibr_parameters[,1] <- Calibrated_parameters

#TEST DELETE AFTERWARDS

Calibr_parameters[3,1] <- Calibr_parameters[3,1]*1.1

source("R calibration\\Main_model_calibration.R") # run the model

output <- f.calibr.output(results_no_screen) # extract the outputs
Predict <- cbind(output$rate_outcomes_m, output$rate_outcomes_f)

# log likelihood instead of sum of squared errors
#m.GOF <- f.GOF.calc(run, m.GOF, Targets, SE, Predict)


c.names <- colnames(Predict)

Predict <- cbind(as.numeric(rownames(Predict)),Predict)
  
colnames(Predict) <- c("Age",c.names)
  
# Extract the calibration targets and calculate the confidence interval
# The incidence is already converted per 1 person alive
Targets <- read.table("R calibration/Calibration_targets.txt", header = TRUE, row.names=1) # Load targets by 5-years
Targets <- f.targets.per.alive.long(Targets, start_age=30, end_age =95, n_by=5) #Linearly extrapolate between the 5-year bands

n_targets <- ncol(Targets) -1 # Exclude Age column
v.target_names <- colnames(Targets)[-1] # Exclude Age column
CI_targets <- f.CI.calc(Targets)

# Calculate SE from the CI
SE_males <- (CI_targets[ ,7:12] - CI_targets[ ,1:6])/3.92
SE_females <- (CI_targets[ ,19:24] - CI_targets[ ,13:18])/3.92
SE <- cbind(SE_males, SE_females)
colnames(SE) <- paste(colnames(Targets)[-1],"_SE", sep="")

# Data sets
Targets <- as.data.frame(Targets)
CI_targets <- as.data.frame(CI_targets)

# Visualize 
# Total incidence rate for males

(target_i.males = ggplot(data = Targets,
                         aes(x = Age,
                             y = incidence.rate.male)) +
    geom_errorbar(aes(ymin = CI_targets$incidence.rate.male_lb, ymax = CI_targets$incidence.rate.male_ub)) +
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

# Stage 3,4 incidence males
(target_i.3.4.males = ggplot(data = Targets,
                             aes(x = Age,
                                 y = inc3.4.male)) +
    geom_errorbar(aes(ymin = CI_targets$inc3.4.male_lb, ymax = CI_targets$inc3.4.male_ub)) +
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
                               y = incidence.rate.fem)) +
    geom_errorbar(aes(ymin = CI_targets$incidence.rate.fem_lb, ymax = CI_targets$incidence.rate.fem_ub)) +
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
(target_i.3.4.females = ggplot(data = Targets,
                               aes(x = Age,
                                   y = inc3.4.fem)) +
    geom_errorbar(aes(ymin = CI_targets$inc3.4.fem_lb, ymax = CI_targets$inc3.4.fem_ub)) +
    geom_point() +
    geom_line(data= as.data.frame(Predict), aes(x = Age, y = BC_stage3.4_female ), color = 'red') +
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
