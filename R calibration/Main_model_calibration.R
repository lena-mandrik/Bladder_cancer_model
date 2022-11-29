# Lena Mandrik 12.07.2022
# This is the main script
###Run this script to run whole model

### run the model:
results_no_screen = foreach::foreach(iterator = 1:n.loops, .options.RNG = optsN) %dorng% {
  gc()
  
  #Progress bar:
  if(f.get_os() != "windows") {
    setTxtProgressBar(pb, iterator) 
  }
  
#Select appropriate parameter set to use
p.set <- ifelse(run_mode =="PSA", i, 1) 

#Set up the model parameters according to current parameter set
set_parameters(p.set)

# Allocate the time to stage at diagnosis for each person in HSE
m.BC.T.to.Stage <- f.stage.assign(m.BC.T.to.Stage) ##!! Needs to be replaced later to avoid loops (apply or map)

# the model baseline population risk is 0 at the model start (age 30 years)

pop <- f.risk.calc(population) #calculates relative and absolute risks of cancer onset

#Set up random number array for each individual
m.Rand <- generate_random()

Simulate_NHD(n.i, n.t, pop)

}

if(f.get_os() == "windows") {
  stopCluster(cluster)
} else {
  close(pb)
  stopCluster(cluster)
}
