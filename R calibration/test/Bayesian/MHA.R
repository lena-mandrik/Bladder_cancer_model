##  This is the code to run calibration either with a Metropolis Hasting algorithm

# Load starting parameters and set up the step
fitted_params <- (read.table("R calibration/Outputs/parameters_fit.txt", header = F, row.names=1)) #exclude the parameter on RR death for no smokers, as it was calculated

fitted_params <-as.matrix(fitted_params[c("P.onset", "P.onset_low.risk", "P.onset_age", "P.onset_sex",
                                "P.sympt.diag_LGBC", "P.sympt.diag_A_HGBC", "P.sympt.diag_B_HGBC", "P.sympt.diag_Age80_HGBC",
                                "shape.t.StI.StII", "shape.t.StII.StIII", "shape.t.StIII.StIV"), 1])

fitted_params <-rbind(fitted_params, 0.0048) # add a value for P.LGtoHGBC that hasn't been calirbated with the random approach

Calibr_parameters <- f.load.start.param(fitted_params, Calibr_parameters)


############################################################################################
######## Metropolis algorithm ################


  chain = array(dim = c(n_samples+1,nrow(Calibr_parameters), 2), dimnames = list(c(1:(n_samples+1)), rownames(Calibr_parameters), c("Param1", "Param2")))
  
  report_prob = matrix(nrow=n_samples+1, 0, ncol=1)
  
  posterior_rep =matrix(nrow=n_samples+1, ncol=2)
  
  # Record the parameters
  chain[1, ,1] = Calibr_parameters[,"Param1"] #Param1
  chain[1, ,2] =  Calibr_parameters[,"Param2"] #Param2
  
  
  for (i in 1:n_samples){
    
    posterior_current = posterior(Calibr_parameters) #Calculate posterior with the initial set and priors
    
    #proposal = proposalfunction(t(chain[i,]), step)
    
    Calibr_parameters <<- f.proposal.param(Calibr_parameters) #sample the parameters 
    
    #Calibr_parameters[,1] <<- proposal #save to the Calibr_parameters in Global env
    posterior_proposal = posterior(Calibr_parameters)
    
    probab = exp(posterior_proposal - posterior_current)
    
    if (runif(1) < probab | runif(1) <0.05){ #accept all parameters with better fit and 5% with worse fit
      
      # Replace all the parameters and record them
      Calibr_parameters <<- f.undate.means(Calibr_parameters, Param_sets)
      chain[(i+1), ,1] <- Calibr_parameters[,"Param1"] #record the sampled parameters
      chain[(i+1), ,2] <- Calibr_parameters[,"Param2"] #record the sampled parameters
        
      report_prob[i+1]=1
      
    }else{
      chain[(i+1), ,1:2] <- chain[i, ,1:2]
    }
    
    posterior_rep[(i+1),1]=posterior_proposal
    posterior_rep[(i+1),2]=posterior_current
    
  
  # If the acceptance is low, decrease the step
  if(i%%50==0){
    if(sum(report_prob[(i-50):i])/50 < 0.1){Calibr_parameters[ ,"step"] <<- Calibr_parameters[ ,"step"]*0.9} #decrease step if low acceptance rate
    
  } #end if
 
    if(i%%10==0){ #Save each 10 loops
      # save intermediary results
      Param1 <- chain[, ,1]
      Param2 <- chain[, ,2]
      write.csv(Param1, file="R calibration\\Outputs\\Param1.csv")
      write.csv(Param2, file="R calibration\\Outputs\\Param2.csv")
      output <- cbind(report_prob,posterior_rep)
      colnames(output) <- c("Accept", "posterior_proposal", "posterior_current")
      write.csv(output, file="R calibration\\Outputs\\output_calibration.csv")
    } #end if
  
  } #end loops

print(chain)