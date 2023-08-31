### Code for processing the results that come out from the simulation ###

#' @details
#' Function for calculating discounted, cumulative, per person results from outcomes matrix 
#' @params
#' results: named list of results
#' @return discounted, cumulative, per person results
#' 
#' For Deterministic analysis
process_DA_results <- function(results, output_file){
  #Sum data from each loop and divide by number of loops 
  ann_results <- Reduce('+', results) / (length(results))
  #Calculate cumulative results
  cum_results <- t(apply(ann_results,1,cumsum))
  #Write data table
  write.table(cum_results, output_file)
  #Return the data table
  cum_results
}

#' 
#' For Probabilistic analysis
process_PSA_results <- function(results, output_file){
  #Convert list into 3D array
  psa_results <- lapply(results, rowSums)
  psa_results <- do.call(cbind, psa_results)
  #Write data table
  write.table(psa_results, output_file)
  #Return the data table
  psa_results
}


#' For testing mode. Analysis of the results
#' 

#list_results <- results_screen_70

f.test.NHD <- function(list_results){

l.Diag <- map(list_results, ~.x$m.Diag)

#Get time from onset to diagnosis for those who get diagnosed during their lifetime

outcomes <- map(l.Diag, f.testing.output)

# extract and unlist the elements
#time to HG diagnosis
l.time.all.diagnosed_HG <- map(outcomes, ~.x$time.all.diagnosed_HG)
time.all.diagnosed_HG <- unlist(l.time.all.diagnosed_HG)
               
l.time.all.diagnosed_HG_sc <- map(outcomes, ~.x$time.all.diagnosed_HG_sc)
time.all.diagnosed_HG_sc <- unlist(l.time.all.diagnosed_HG_sc)             

l.time.all.diagnosed_HG_sy <- map(outcomes, ~.x$time.all.diagnosed_HG_sy)
time.all.diagnosed_HG_sy <- unlist(l.time.all.diagnosed_HG_sy)

#time to LG diagnosis
l.time.all.diagnosed_LG <- map(outcomes, ~.x$time.all.diagnosed_LG)
time.all.diagnosed_LG <- unlist(l.time.all.diagnosed_LG)

l.time.all.diagnosed_LG_sc <- map(outcomes, ~.x$time.all.diagnosed_LG_sc)
time.all.diagnosed_LG_sc <- unlist(l.time.all.diagnosed_LG_sc)             

l.time.all.diagnosed_LG_sy <- map(outcomes, ~.x$time.all.diagnosed_LG_sy)
time.all.diagnosed_LG_sy <- unlist(l.time.all.diagnosed_LG_sy) 

# Time to death
l.time.to.death <- map(outcomes, ~.x$time.to.death)
time.to.death <- unlist(l.time.to.death)

l.time.to.death_sy <- map(outcomes, ~.x$time.to.death_sy)
time.to.death_sy <- unlist(l.time.to.death_sy)

l.time.to.death_sc <- map(outcomes, ~.x$time.to.death_sc)
time.to.death_sc <- unlist(l.time.to.death_sc)

l.LG.progr.HG <- map(outcomes, ~.x$LG.progr.HG)
LG.progr.HG <- sapply(l.LG.progr.HG, mean)
LG.progr.HG[is.nan(LG.progr.HG)] <-0

mean.LG.progr.HG <- mean(LG.progr.HG)

}


#####################################
# Create a histogram plot
#plot_time_to_event <- function(time_to_event, main, xlab, ylab) {
  # Create a histogram plot
 # hist(time_to_event, 
  #     xlab = xlab,
   #    ylab = ylab,
   #    main = main)
  
  # Calculate mean and standard deviation
 # mean_value <- mean(time_to_event)
 # sd_value <- sd(time_to_event)
  
  # Add mean and standard deviation to the plot
  #abline(v = mean_value, col = "red", lwd = 2)
 # text(mean_value, max(hist(time_to_event)$counts), "Mean", col = "red")
 # abline(v = mean_value - sd_value, col = "blue")
 # abline(v = mean_value + sd_value, col = "blue")
 # text(mean_value - sd_value, min(hist(time_to_event)$counts), "Mean - SD", col = "blue")
 # text(mean_value + sd_value, max(hist(time_to_event)$counts), "Mean + SD", col = "blue")
  
  
#}

#plot.mort = plot_time_to_event(time.to.death, main = "Distribution of Time to Death", xlab = "Time to Death (Years)", ylab = "Frequency")


#hist(time.all.diagnosed_LG, 
   #  xlab = "Time to diagnosis, LR (Years)",  # Set x-axis label
   #  ylab = "Frequency",              # Set y-axis label
   #  main = "Distribution of Time to Diagnoses")  # Set plot title

# Calculate mean and standard deviation
#mean_value <- mean(time.all.diagnosed_LG)
#sd_value <- sd(time.all.diagnosed_LG)

# Add mean and standard deviation to the plot
#abline(v = mean_value, col = "red", lwd = 2)    # Add vertical line for mean
#text(mean_value, 4, "Mean", col = "red")         # Add label for mean
#abline(v = mean_value - sd_value, col = "blue")  # Add vertical line for mean - standard deviation
#abline(v = mean_value + sd_value, col = "blue")  # Add vertical line for mean + standard deviation
#text(mean_value - sd_value, 6, "Mean - SD", col = "blue")  # Add label for mean - standard deviation
#text(mean_value + sd_value, 6, "Mean + SD", col = "blue")  #

#t_test <- t.test(time.all.diagnosed_LG)
#ci_lower <- t_test$conf.int[1]
#ci_upper <- t_test$conf.int[2]

# Print the results
#cat("95% Confidence Interval:", ci_lower, "to", ci_upper)
##############################################

#' Function returning time to diagnosis, time to death from onset, and proportion of LG progressing to HG
#' 
f.testing.output <- function(m.Diag){
  
  # all diagnosed pop
  all.diagnosed_HG <- m.Diag[m.Diag[ ,"HG_stage_diag"]>0,]
  all.diagnosed_LG <- m.Diag[m.Diag[ ,"LG_diag"]>0,]
  
  #Get time from onset to diagnosis for those who get diagnosed during their lifetime
  time.all.diagnosed_HG <- all.diagnosed_HG[ , "HG_age_diag"] - all.diagnosed_HG[ , "HG_age_onset"]
  time.all.diagnosed_LG <- all.diagnosed_LG[ , "LG_age_diag"] - all.diagnosed_LG[ , "LG_age_onset"]
  
  # time to death for those who died from BC
  all.died.BC <- all.diagnosed_HG[all.diagnosed_HG[,"age_BC_death"]>0,]
  time.to.death <- all.died.BC[ , "age_BC_death"] - all.died.BC[ , "HG_age_onset"]
  
  # for sympt diagnosed
    # all diagnosed pop
    all.diagnosed_HG_sy <- m.Diag[m.Diag[ ,"HG_stage_diag"]>0 & m.Diag[ ,"HG_sympt_diag"]==1,]
    all.diagnosed_LG_sy <- m.Diag[m.Diag[ ,"LG_diag"]>0& m.Diag[ ,"LG_sympt_diag"]==1,]
  
    #Get time from onset to diagnosis for those who get diagnosed during their lifetime
    time.all.diagnosed_HG_sy  <- all.diagnosed_HG_sy[ , "HG_age_diag"] - all.diagnosed_HG_sy[ , "HG_age_onset"]
    time.all.diagnosed_LG_sy  <- all.diagnosed_LG_sy[ , "LG_age_diag"] - all.diagnosed_LG_sy[ , "LG_age_onset"]
  
    # time to death for those who died from BC
    all.died.BC_sy <- all.diagnosed_HG_sy[all.diagnosed_HG_sy[,"age_BC_death"]>0,]
    time.to.death_sy <- all.died.BC_sy[ , "age_BC_death"] - all.died.BC_sy[ , "HG_age_onset"]
    
  # for screen diagnosed
    # all diagnosed pop
    all.diagnosed_HG_sc <- m.Diag[m.Diag[ ,"HG_stage_diag"]>0 & m.Diag[ ,"HG_screen_diag"]==1,]
    all.diagnosed_LG_sc <- m.Diag[m.Diag[ ,"LG_diag"]>0& m.Diag[ ,"LG_screen_diag"]==1,]
    
    if(is.vector(all.diagnosed_HG_sc)==TRUE){time.all.diagnosed_HG_sc  <- all.diagnosed_HG_sc["HG_age_diag"] - all.diagnosed_HG_sc["HG_age_onset"]} else{
      #Get time from onset to diagnosis for those who get diagnosed during their lifetime
      time.all.diagnosed_HG_sc  <- all.diagnosed_HG_sc[ , "HG_age_diag"] - all.diagnosed_HG_sc[ , "HG_age_onset"]}
    
    if(is.vector(all.diagnosed_LG_sc)==TRUE){time.all.diagnosed_LG_sc  <- all.diagnosed_LG_sc["LG_age_diag"] - all.diagnosed_LG_sc["LG_age_onset"]} else{
      time.all.diagnosed_LG_sc  <- all.diagnosed_LG_sc[ , "LG_age_diag"] - all.diagnosed_LG_sc[ , "LG_age_onset"]
    }
    
  
    
    # time to death for those who died from BC
    #all.died.BC_sc <- all.diagnosed_HG_sc[all.diagnosed_HG_sc[,"age_BC_death"]>0,]
    #if(is.vector(all.died.BC_sc)==TRUE){time.to.death_sc <- all.died.BC_sc["age_BC_death"] - all.died.BC_sc["HG_age_onset"]} else{
    #  time.to.death_sc <- all.died.BC_sc[ , "age_BC_death"] - all.died.BC_sc[ ,"HG_age_onset"]
      
   # }

    # Number of people with LG who progressed to HG and the total LG 
    LG.to.HG <- nrow(m.Diag[m.Diag[,"LG_state"]==1 & m.Diag[,"HG_state"]==1, ])
    All.LG <- nrow(m.Diag[m.Diag[,"LG_state"]==1, ])
    LG.progr.HG <- LG.to.HG/All.LG
    
    l.output = list(time.all.diagnosed_HG=time.all.diagnosed_HG,
                    time.all.diagnosed_HG_sy=time.all.diagnosed_HG_sy,
                    time.all.diagnosed_HG_sc=time.all.diagnosed_HG_sc,
                    time.all.diagnosed_LG=time.all.diagnosed_LG,
                    time.all.diagnosed_LG_sy=time.all.diagnosed_LG_sy,
                    time.all.diagnosed_LG_sc=time.all.diagnosed_LG_sc,
                    time.to.death=time.to.death,
                    time.to.death_sy=time.to.death_sy,
                    #time.to.death_sc=time.to.death_sc,
                    LG.progr.HG=LG.progr.HG
                    )
    return(l.output)
}

