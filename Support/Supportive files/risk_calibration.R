# The calibration of the risks for smokers and past smokers

f.risk.calibration <- function(results_no_screen){
  
  l.Diag <- map(results_no_screen, ~.x$m.Diag)
  #TR_m <- map(results_no_screen, ~.x$TR_m)
  #TR_f <- map(results_no_screen, ~.x$TR_f)
  
  #Set up the model parameters according to current parameter set
  f.set_parameters(1)
  
  pop <- f.risk.calc(population) 
  
  incidence.smoke <- map(l.Diag, f.diag.outcomes)
  outcomes_smoke <- (Reduce('+', incidence.smoke) / (length(incidence.smoke)))
  
  pop_smoke <- nrow(subset(pop, pop[ , "current_smoke"] ==1))
  pop_past.smoke <- nrow(subset(pop, pop[ , "past_smoke"] ==1))
  pop_no.smoke <- nrow(subset(pop, pop[ , "no_smoke"] ==1))
  
  pop_smoke <- c(pop_past.smoke, pop_smoke, pop_no.smoke)
  
  incidence_by_smoke <- outcomes_smoke/pop_smoke
  
  RR_smoke <- incidence_by_smoke[2]/incidence_by_smoke[3]
  RR_past.smoke <- incidence_by_smoke[1]/incidence_by_smoke[3]
  
  risks <- cbind(RR_smoke, RR_past.smoke)
  
  risks
  
}


f.diag.outcomes <- function(m.Diag){
  
  pop1 <- pop[ , c("current_smoke", "past_smoke", "no_smoke")]
  m.Diag <- cbind(m.Diag, pop1)
  
  m.diagnosed <- subset(m.Diag, m.Diag[ , "BC_diag"] ==1, select = c("BC_diag", "current_smoke", "past_smoke", "no_smoke"), drop = FALSE)
  
  smoke <- nrow(subset(m.diagnosed, m.diagnosed[ , "current_smoke"] ==1))
  past_smoke <- nrow(subset(m.diagnosed, m.diagnosed[ , "past_smoke"] ==1))
  no_smoke <- nrow(subset(m.diagnosed, m.diagnosed[ , "no_smoke"] ==1))
  
  n_events <- c(smoke, past_smoke, no_smoke)
  
  n_events
}