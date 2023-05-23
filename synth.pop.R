# Function not called in the model but used to develop the modelling input
#' @details
#' This function to calculate probability of population be in each disease state based on their age at the HSE
#' @params
#' 
#' 
#' 
#' @return matrix with probabilities for each person
#' 
library("plyr") #required library

f.states.prob <- function(list_results){
  
  l.Diag <- map(list_results, "m.Diag") #matrix to extract cancer incidence
  l.M <- map(list_results, "m.M_8s") #matrix to extract probability of states
  
  m.M <- do.call("rbind",l.M)
  m.Diag <- do.call("rbind",l.Diag)
  
  
  m.Diag <- (cbind( m.M[ ,1:3],m.Diag))
  colnames(m.M) <- c("PID", "age", "cycle", 1:71)
  
  names_m.Diag <- colnames(m.Diag)[-c(1,2)]
  colnames(m.Diag) <- c("PID", "age", names_m.Diag)
  
  m.Diag <- m.Diag[ ,c("PID", "age", "cycle", "LG_age_diag", "HG_age_diag", "HG_stage_diag", "HG_diag" , "LG_diag")]
  m.Diag[ , "LG_age_diag"] <- replace( m.Diag[ , "LG_age_diag"],  m.Diag[ , "LG_age_diag"]==0, 100)
  m.Diag[ , "HG_age_diag"] <- replace( m.Diag[ , "HG_age_diag"],  m.Diag[ , "HG_age_diag"]==0, 100)
  
  
  # probability for age x
  m.pop.prob <- matrix(0,nrow=n.i, ncol =n.s_long)
  colnames(m.pop.prob) <- 1:8
  m.BC <- matrix(0,nrow=n.i, ncol =5)
  colnames(m.BC) <- 1:5
  
  m.diagnosed <- matrix(0,nrow=n.i, ncol =5)
    
  for(i in 1:n.i){
    cycle <- m.M[i,"cycle"]+3 #number of column to look for
    
   # if(m.M[i, "age"]>age1 & m.M[i, "age"] <age2){ # condition to do this only for pop of screening age
     # if(m.M[i, "age"]>30){

      for(n in 1:n.s_long){
        count_events <-  as.matrix(count(m.M[m.M[,"PID"]==i,cycle]))
        condition <- count_events[count_events[,1]== n,1]
        m.pop.prob[i,n] <- replace(m.pop.prob[i,n], n==condition, count_events[count_events[,1]==n,2])
        
           }
     #else{m.pop.prob[i,1]<-1} #for if condition, assume NE at age 30y
    m.BC[i,1] <- sum(((m.Diag[m.M[,"PID"]==i & m.Diag[,"LG_age_diag"] < m.M[ ,"age"] & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "LG_age_diag"] - m.Diag[m.Diag[,"PID"]==i & m.Diag[,"LG_age_diag"]<m.M[ ,"age"] & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "age"]) <0)*1)
    m.BC[i,2] <- sum(((m.Diag[m.M[,"PID"]==i & m.Diag[,"HG_stage_diag"]==1 & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "HG_age_diag"] - m.Diag[m.Diag[,"PID"]==i & m.Diag[,"HG_stage_diag"]==1 & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "age"]) <0)*1)
    m.BC[i,3] <- sum(((m.Diag[m.M[,"PID"]==i & m.Diag[,"HG_stage_diag"]==2 & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "HG_age_diag"] - m.Diag[m.Diag[,"PID"]==i & m.Diag[,"HG_stage_diag"]==2 & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "age"]) <0)*1)
    m.BC[i,4] <- sum(((m.Diag[m.M[,"PID"]==i & m.Diag[,"HG_stage_diag"]==3 & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "HG_age_diag"] - m.Diag[m.Diag[,"PID"]==i & m.Diag[,"HG_stage_diag"]==3 & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "age"]) <0)*1)
    m.BC[i,5] <- sum(((m.Diag[m.M[,"PID"]==i & m.Diag[,"HG_stage_diag"]==4 & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "HG_age_diag"] - m.Diag[m.Diag[,"PID"]==i & m.Diag[,"HG_stage_diag"]==4 & m.M[, cycle] !=4 & m.M[, cycle] !=8,  "age"]) <0)*1)
    
  }
  
  length <- rowSums(m.pop.prob[, 1:8])
  
  pop.prob <- m.pop.prob[,1:8]/length
  m.BC.diag <- m.BC/length
  
  # Calculate a probability to be a diagnosed BC patient by HSE age
  symptomatic <- pop.prob[, c(2,3,5:7)]*m.BC.diag
  all.symptomatic <- rowSums(symptomatic)
    
  # Calculate a probability to be in a non-diagnosed BC state
  assymptomatic <-pop.prob[, c(2,3,5:7)] - symptomatic
  
  # final states prob
  states.prob <- pop.prob
  states.prob[, c(2,3,5:7)] <- assymptomatic
  states.prob <- cbind(states.prob, all.symptomatic)
  
  write.table(states.prob,"probabilities_of_states.txt", sep = "\t")
  
  return(states.prob)
}



pop_prob <- f.states.prob(results_no_screen)


# adjust the probabilities to ignore the probability to die (since everyone in the HSE is alive)

states_pop <- as.matrix(read.table("Data\\probabilities_of_states.txt", header = TRUE))
colnames(states_pop) <- c(states_long, "diagnosed")

no_BC <- rowSums(cbind(states_pop[,"NoBC"],states_pop[,"DeathOC"]))
diagnosed <- rowSums(cbind(states_pop[,"DeathBC"],states_pop[,"diagnosed"]))

# replace with the probabilities combined for those who died from BC and other causes

states_pop[ ,"NoBC"] = no_BC
states_pop[ ,"DeathOC"] = 0
states_pop[ ,"DeathBC"] = 0
states_pop[ ,"diagnosed"] = diagnosed
#states_pop[is.na(states_pop[ ,"NoBC"]),"NoBC"] <-1
#states_pop[is.na(states_pop[ ,"BC_LG"]),2:9] <-0

write.table(states_pop, "Data\\states_pop.txt", col.names=T)
