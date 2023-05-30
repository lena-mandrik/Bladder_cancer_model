
#check the code
#check whether the matrices m.M and m.M_8 are identical

m.check = matrix(0, nrow=ncol(m.M), ncol=8)
colnames(m.check) <- c("NC_m.M", "NC_m.M8", "LG_m.M", "LG_m.M8", "HG_m.M", "HG_m.M8", "death_m.M", "death_m.M8")

for (i in 1:ncol(m.M)){
  
  #no cancer
  m.check[i,1] =length(which(m.M[,i]==1))
  m.check[i,2] =length(which(m.M_8s[,i]==1))
  
  #LG
  m.check[i,3] =length(which(m.M[,i]==2))
  m.check[i,4] =length(which(m.M_8s[,i]==2))
  
  #HG
  m.check[i,5] =length(which(m.M[,i]==3))
  m.check[i,6] =length(which(m.M_8s[,i]==3| m.M_8s[,i]==5 |m.M_8s[,i]==6 |m.M_8s[,i]==7))
  
  #death
  
  m.check[i,7] =length(which(m.M[,i]==4))
  m.check[i,8] =length(which(m.M_8s[,i]==4| m.M_8s[,i]==8))
}

#check the pop  wtih HG onset has HG state

cancer_m.Diag =which(m.Diag[ ,"HG_state"]==1 |m.Diag[ ,"LG_state"]==1)
HG_onset = which(m.Diag[ ,"HG_yr_onset"]>0)
HG_state = which(m.Diag[ ,"HG_state"]>0)

identical(HG_onset, HG_state)

# check mort from BC is reflected 
m.Diag =cbind(m.Diag,pop[ ,"PID"])
m.Diag_c = m.Diag[cancer_m.Diag, ]

p15_a=83
p214_a=67
p448_a=71
p500_a=93
p611_a=94
p668_a=84

pop_cd=pop[c(3,41,47),]

######################################
#for validation
m.Diag_ns <- results_no_screen[[1]]
m.Diag_s <- results_screen_70[[1]]

#Diagnosed
LG_sc = length(which(m.Diag_ns[ ,"LG_screen_diag"]==1))
HG_sc = length(which(m.Diag_ns[ ,"HG_screen_diag"]==1))
LG_sym = length(which(m.Diag_ns[ ,"LG_sympt_diag"]==1))
HG_sym = length(which(m.Diag_ns[ ,"HG_sympt_diag"]==1))
NS_total = LG_sc+ LG_sym + HG_sc  + HG_sym

sLG_sc = length(which(m.Diag_s[ ,"LG_screen_diag"]==1))
sHG_sc = length(which(m.Diag_s[ ,"HG_screen_diag"]==1))
sLG_sym = length(which(m.Diag_s[ ,"LG_sympt_diag"]==1))
sHG_sym = length(which(m.Diag_s[ ,"HG_sympt_diag"]==1))
S_total = sLG_sc+ sLG_sym + sHG_sc  + sHG_sym

# onset
(LG_state_NS = length(which(m.Diag_ns[ ,"LG_state"]==1)))
(HG_state_NS = length(which(m.Diag_ns[ ,"HG_state"]==1)))

(LG_state_S = length(which(m.Diag_s[ ,"LG_state"]==1)))
(HG_state_S = length(which(m.Diag_s[ ,"HG_state"]==1)))

