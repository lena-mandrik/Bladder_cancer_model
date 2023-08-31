###############################################
RR_smoke =2.5 #1.3900000
RR_past.smoke =1.20000

m.risks =matrix(0, nrow=100, ncol=4)
colnames(m.risks) =c("RR_smoke",  "RR_past.smoke", "RR_smoke_predict",  "RR_past.smoke_predict")

for(x in 1:nrow(m.risks)){
  
  m.risks[x,1:2] =c(RR_smoke, RR_past.smoke)
  
  source("Main_model.R")
  
  risk = f.risk.calibration(results_no_screen, pop)
  m.risks[x,"RR_smoke_predict"] = risk[1]
  m.risks[x, "RR_past.smoke_predict"] =risk[2]
  
if(abs((risk[1] -1.390000e+00)/1.390000e+00) >0.001){RR_smoke=RR_smoke*1.01} 
  if(abs((risk[1] -1.2)/1.2) >0.001){RR_past.smoke=RR_past.smoke*1.01} 
  
}


write.csv(m.risks, file="KC_risk_calibr.csv")