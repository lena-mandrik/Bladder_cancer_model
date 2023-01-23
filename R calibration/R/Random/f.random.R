##  This is the code to run calibration either with a random approach 



# Search dimensions. Use for random search algorithm only
lower_boud <- c(0.000001, 0.1, 1.001, 1.001, 0.005, 0.001, 1.01, 0.85, -0.3, rep(1.005,3)) #for symtomatic presentation requires other limits - A <15%, B =5-40%, C=10-60%, D-30-100%
upper_boud <- c(0.001, 0.6, 1.2, 1.9, 0.1, 0.2, 1.9, 1, -0.01, rep(1.5,3))

# Sample using Latin Hypercube
sample_LH <- randomLHS(n_samples, n_params)

# Re-scale to min/max
m.sample.params <- matrix(nrow=n_samples, ncol = n_params)
colnames(m.sample.params) <- v.param_names

for(i in 1:n_params){
  m.sample.params[,i] <- qunif(sample_LH[,i],
                               min = lower_boud[i],
                               max = upper_boud[i])
}
write.csv(m.sample.params, file="R calibration\\ParametersInput.csv")




# Run the calibration loops

for(run in 1:n_samples){
  
  Calibr_parameters <- as.matrix(m.sample.params[run,], ncol = 1) # get the inputs for calibrating parameters
  
  source("R calibration\\Main_model_calibration.R") # run the model
  
  output <- f.calibr.output(results_no_screen) # extract the outputs
  
  Predict <- cbind(output$rate_outcomes_m, output$rate_outcomes_f)
  
  # log likelihood instead of sum of squared errors
  
  m.GOF <- f.GOF.calc(run, m.GOF, Targets, SE, Predict)
  
  
  if(run%%2==0){
    cat('\r', paste(round(run/n_samples*100), "% done", sep = " "))
    write.csv(m.GOF, file="R calibration\\Outputs\\m.GOF.csv")}
} 
output <- cbind(m.sample.params, m.GOF)



