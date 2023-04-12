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


