##  This is the code to run calibration either with a Metropolis Hasting algorithm


start <- Sys.time()

MHA.step =0.3


Calibr_parameters["P.onset",1]=5.5E-06 
Calibr_parameters["P.onset_low.risk",1]= 0.6376
Calibr_parameters["P.onset_age",1]=1.138
Calibr_parameters["RR.onset_sex",1]=3.64 

Calibr_parameters["P.sympt.diag_LGBC",1]=0.102
Calibr_parameters["P.sympt.diag_St1",1]=0.15
Calibr_parameters["P.sympt.diag_St2",1]=0.25
Calibr_parameters["P.sympt.diag_St3",1]=0.33 
Calibr_parameters["P.sympt.diag_St4",1]=0.34 
Calibr_parameters["P.sympt.diag_Age",1]= 0.91
Calibr_parameters["shape.t.StI.StII",1]=1
Calibr_parameters["shape.t.StII.StIII",1]= 1
Calibr_parameters["shape.t.StIII.StIV",1]= 1
Calibr_parameters["P.LGtoHGBC",1]= 2.555e-03
Calibr_parameters["P.ungiag.dead",1]=-0.0051 
Calibr_parameters["Mean.t.StI.StII",1]= 6.5#5.2 #4.8
Calibr_parameters["Mean.t.StII.StIII",1]=4.2
Calibr_parameters["Mean.t.StIII.StIV",1]=4

############################################################################################
######## Metropolis algorithm ################


  chain = matrix(0, nrow= n_samples+1, ncol =nrow(Calibr_parameters)+3) #array(dim = c(n_samples+1,nrow(Calibr_parameters), 2), dimnames = list(c(1:(n_samples+1)), rownames(Calibr_parameters), c("Param1", "Param2")))
  colnames(chain) = c(rownames(Calibr_parameters), "accept", "posterior_current", "posterior_proposal")
                      

   chain[1, 1:nrow(Calibr_parameters)]=  Calibr_parameters[,"Mean"]
   
   posterior_current = posterior(Calibr_parameters) #Calculate posterior with the initial set and priors
  
  for (i in 1:n_samples){
    
    Calibr_parameters <<- f.proposal.param(Calibr_parameters, MHA.step) #sample the parameters 
    
    posterior_proposal = posterior(Calibr_parameters)
    
    probab = exp(posterior_proposal - posterior_current)
    
    if (runif(1) < probab |(i>300 & (runif(1) < 0.05 & (posterior_proposal - posterior_current)/posterior_current < 0.5))){ #accept all parameters with better fit and 1% with worse fit
      
      chain[(i+1), 1:nrow(Calibr_parameters)] =  Calibr_parameters[,"Mean"]
      
      chain[(i+1), "accept"]=1
      
      posterior_current =posterior_proposal

    }else{
      chain[(i+1), 1:nrow(Calibr_parameters)]= chain[(i), 1:nrow(Calibr_parameters)]
      Calibr_parameters[,1] <- chain[(i+1), 1:nrow(Calibr_parameters)] #update the parameters depending whether they were accepted or not
          }
    
    chain[(i+1), "posterior_proposal"] = posterior_proposal
    chain[(i+1), "posterior_current"] = posterior_current
    
  
  # If the acceptance is low, decrease the step
  if(i%%100==0){
    if(sum(chain[(i-100):i, "accept"])/100 < 0.1){MHA.step <- MHA.step*0.9} #decrease step if low acceptance rate
    
  } #end if
 
    if(i%%10==0){ #Save each 10 loops

      write.csv(chain, file="R calibration\\Outputs\\MHA_chain.csv")
      
      cat('\r', paste(round(i/n_samples*100), "% done", sep = " "))
    } #end if
  
  } #end loops

print(chain)

print(Sys.time() -start)