### CODE TO PROCESS HSE2014 DATA TXT FILE ###
### Chloe Thomas May 2017 ##
### Lena Mandrik January 2019#

# Load up dataset
Dataset<-as.matrix(read.table("X:\\ScHARR\\PR_Opt_Bowel_Cancer_Screening\\General\\PHASE 2 Optimising bowel screening\\Model\\Current Model Data\\Population HSE2014\\UKDA-7919-tab\\tab\\hse2014ai.txt",header=TRUE))
newdata<-matrix(0,ncol=20,nrow=length(Dataset[,1]))

# Select columns to keep (column names from data dictionary)
newdata[,1]<-Dataset[,'pserial'] # HSE serial number, Numeric
newdata[,2]<-Dataset[,'Sex'] # Sex, Numeric 1 = male; 2 = female
newdata[,3]<-Dataset[,'Age90'] # Age, Numeric
newdata[,4]<-Dataset[,'qimd'] # Deprivation score, Numeric 1 = least deprived; 5 = most deprived
newdata[,5]<-Dataset[,'origin2'] # Ethnicity. Numeric 1 = white; 2 = black; 3 = Asian; 4 = mixed; 5 = other
newdata[,6]<-Dataset[,'BMIval'] # BMI. Numeric
newdata[,7]<-Dataset[,'cigsta3'] # Cigarette Smoking Status. Numeric 1= current; 2 = Ex-regular smoker; 3 = Never regular smoker
newdata[,9]<-Dataset[,'dnevr'] # Whether always non-drinker. Numeric 1 = always non-drinker; 2 = used to drink but stopped.
newdata[,10]<-Dataset[,'totalwu'] # Total units of alcohol in week. Numeric
newdata[,11]<-Dataset[,'TotmModWk'] # Time spent on moderate PA in the last week (minutes)
newdata[,12]<-Dataset[,'TotmVigWk'] # Time spent on vigorous PA in the last week (minutes)
newdata[,14]<-Dataset[,'Mobility'] # EQ-5D response to mobility
newdata[,15]<-Dataset[,'Selfcare'] # EQ-5D response to selfcare
newdata[,16]<-Dataset[,'UsualAct'] # EQ-5D response to Usual activity
newdata[,17]<-Dataset[,'Pain'] # EQ-5D response to pain
newdata[,18]<-Dataset[,'Anxiety'] # EQ-5D response to anxiety
newdata[,20]<-Dataset[,'wt_int'] # Weight for analysis of core interview sample

colnames(newdata)<-c("PID","sex","age_0","imd","ethnic","bmi","smoke","smo","alc_hist","units", "ModPA", "VigPA","PA","Mobility","Selfcare","UsualAct","Pain","Anxiety","EQ5D","weighting")

##Exclude data from individuals aged under 30
newdata <- subset(newdata,newdata[,'age_0']>=30)
length(newdata[,1])

##Summary statistics
table(newdata[,'sex'])
summary(newdata[,'age_0'][newdata[,'age_0']>0])
sd(newdata[,'age_0'][newdata[,'age_0']>0])
table(newdata[,'age_0'])
table(newdata[,'ethnic'])
table(newdata[,'imd'])
table(newdata[,'bmi'])[1:10]
summary(newdata[,'bmi'][newdata[,'bmi']>0])
sd(newdata[,'bmi'][newdata[,'bmi']>0])
table(newdata[,'smoke'])
table(newdata[,'alc_hist'])
table(newdata[,'units'])[1:10]
summary(newdata[,'units'][newdata[,'units']>0])
sd(newdata[,'units'][newdata[,'units']>0])
table(newdata[,'ModPA'])[1:10]
summary(newdata[,'ModPA'][newdata[,'ModPA']>0])
sd(newdata[,'ModPA'][newdata[,'ModPA']>0])
table(newdata[,'VigPA'])[1:10]
summary(newdata[,'VigPA'][newdata[,'VigPA']>0])
sd(newdata[,'VigPA'][newdata[,'VigPA']>0])
table(newdata[,'Mobility'])[1:5]
table(newdata[,'Selfcare'])[1:5]
table(newdata[,'UsualAct'])[1:5]
table(newdata[,'Pain'])[1:5]
table(newdata[,'Anxiety'])[1:5]

#Missing data analysis
##replace missing data with NA for later imputation
newdata[,'Mobility']<-replace(newdata[,'Mobility'],newdata[,'Mobility']==-1,NA)    
newdata[,'Mobility']<-replace(newdata[,'Mobility'],newdata[,'Mobility']==-9,NA)
newdata[,'Selfcare']<-replace(newdata[,'Selfcare'],newdata[,'Selfcare']==-1,NA)    
newdata[,'Selfcare']<-replace(newdata[,'Selfcare'],newdata[,'Selfcare']==-9,NA)
newdata[,'UsualAct']<-replace(newdata[,'UsualAct'],newdata[,'UsualAct']==-1,NA)    
newdata[,'UsualAct']<-replace(newdata[,'UsualAct'],newdata[,'UsualAct']==-9,NA)
newdata[,'Pain']<-replace(newdata[,'Pain'],newdata[,'Pain']==-1,NA)    
newdata[,'Pain']<-replace(newdata[,'Pain'],newdata[,'Pain']==-9,NA)
newdata[,'Anxiety']<-replace(newdata[,'Anxiety'],newdata[,'Anxiety']==-1,NA)    
newdata[,'Anxiety']<-replace(newdata[,'Anxiety'],newdata[,'Anxiety']==-9,NA)

##Recode data and replace missing data with assumed data where appropriate
###Recode gender data so that 1 = male and 0 = female  
newdata[,'sex']<-replace(newdata[,'sex'],newdata[,'sex']==2,0)

###ethnicity: assume that all those with missing data (n = 27) are white.
newdata[,'ethnic']<-replace(newdata[,'ethnic'],newdata[,'ethnic']==-9|newdata[,'ethnic']==-8,1)

###smoking: assume that all those with missing data ( n = 20) were never regular smokers
newdata[,'smoke']<-replace(newdata[,'smoke'],newdata[,'smoke']==-9|newdata[,'smoke']==-8,3)
###Set up Smo colomn to separate current smokers (1) from current non-smokers (0)
newdata[,'smo']<-replace(newdata[,'smo'],newdata[,'smoke']==1,1)
table(newdata[,'smo'])
###Modify smoke so that it represents past smokers
newdata[,"smoke"]<-replace(newdata[,'smoke'],newdata[,'smoke']==1|newdata[,'smoke']==3,0)
newdata[,"smoke"]<-replace(newdata[,'smoke'],newdata[,'smoke']==2,1)

###alcohol units in last week: assume that all those with missing data (92) drank average, unless know they are non drinkers 
newdata[,'units']<-replace(newdata[,'units'],(newdata[,'units']==-9|newdata[,'units']==-8) & newdata[,'alc_hist']==-1, mean(newdata[,'units'][newdata[,'units']>0]))
newdata[,'units']<-replace(newdata[,'units'],(newdata[,'units']==-9|newdata[,'units']==-8) & newdata[,'alc_hist']>=1, 0)
table(newdata[,'units'])[1:10]
summary(newdata[,'units'])

##physical activity: assume all those with missing data have no vigorous or moderate PA
newdata[,'VigPA']<-replace(newdata[,'VigPA'],newdata[,'VigPA']==-9 | newdata[,'VigPA']==-1, 0)
newdata[,'ModPA']<-replace(newdata[,'ModPA'],newdata[,'ModPA']==-9 | newdata[,'ModPA']==-1, 0)
#Make a new variable called PA, which corresponds to PA METS
#1 minute of moderate PA corresponds to 4 METS
#1 minute of vigorous PA corresponds to 8 METS
newdata[,'PA']<-newdata[,'ModPA'] * 4 + newdata[,'VigPA'] * 8
table(newdata[,'PA'])[1:10]
summary(newdata[,'PA'])

### EQ5D score
EQ5D<-matrix(0,ncol=12,nrow=length(newdata[,1]))
EQ5D[,1]<-replace(EQ5D[,1],newdata[,'Mobility']==2 |newdata[,'Mobility']==3 |newdata[,'Selfcare']==2 |newdata[,'Selfcare']==3 |newdata[,'Anxiety']==2 |newdata[,'Anxiety']==3 |newdata[,'UsualAct']==2 |newdata[,'UsualAct']==3 |newdata[,'Pain']==2 | newdata[,'Pain']==3,0.081)
EQ5D[,2]<-replace(EQ5D[,2],newdata[,'Mobility']==3 |newdata[,'Selfcare']==3 |newdata[,'Anxiety']==3 |newdata[,'UsualAct']==3 | newdata[,'Pain']==3,0.269)
EQ5D[,3]<-replace(EQ5D[,3],newdata[,'Mobility']==2,0.069)
EQ5D[,4]<-replace(EQ5D[,4],newdata[,'Mobility']==3,0.314)
EQ5D[,5]<-replace(EQ5D[,5],newdata[,'Selfcare']==2,0.104)
EQ5D[,6]<-replace(EQ5D[,6],newdata[,'Selfcare']==3,0.214)
EQ5D[,7]<-replace(EQ5D[,7],newdata[,'UsualAct']==2,0.036)
EQ5D[,8]<-replace(EQ5D[,8],newdata[,'UsualAct']==3,0.094)
EQ5D[,9]<-replace(EQ5D[,9],newdata[,'Pain']==2,0.123)
EQ5D[,10]<-replace(EQ5D[,10],newdata[,'Pain']==3,0.386)
EQ5D[,11]<-replace(EQ5D[,11],newdata[,'Anxiety']==2,0.071)
EQ5D[,12]<-replace(EQ5D[,12],newdata[,'Anxiety']==3,0.236)
newdata[,'EQ5D'] <- 1-rowSums(EQ5D)
newdata[,'EQ5D'] <- replace(newdata[,'EQ5D'],is.na(newdata[,'Mobility']) | is.na(newdata[,'Selfcare']) | is.na(newdata[,'UsualAct']) | is.na(newdata[,'Pain']) | is.na(newdata[,'Anxiety']),-1)

###Code for summarising EQ5D by age, sex and imd, and imputing missing values
# Summarise EQ5D data in individuals that have data
newdata1 <- subset(newdata,newdata[,'EQ5D'] != (-1))
EQ5D_summary <- aggregate(EQ5D ~ sex + imd + age_0, data = newdata1, FUN = mean)
# Linear model for calculating EQ5D - alter values of coefficients below if different
a <- lm(EQ5D ~ sex + imd + age_0, data = EQ5D_summary)
summary(a)
#plot(a)
Intercept <- 1.1041305
sex <- 0.0335711
imd <- (-0.0282924)
age <- (-0.0036214)

# Calculate missing values EQ-5D
imputed_EQ5D <- Intercept + newdata[,'sex']*sex + newdata[,'age_0']*age + newdata[,'imd']*imd
no_EQ5D <- newdata[,'EQ5D']==(-1)
imputed_EQ5D <- imputed_EQ5D * no_EQ5D
newdata[,'EQ5D'] <- replace(newdata[,'EQ5D'], newdata[,'EQ5D']==(-1), 0)
newdata[,'EQ5D'] <- newdata[,'EQ5D'] + imputed_EQ5D
newdata[,'EQ5D'] <- replace(newdata[,'EQ5D'], newdata[,'EQ5D'] >1, 1)
summary(newdata[,"EQ5D"])
sd(newdata[,"EQ5D"])
New_EQ5D_summary <- aggregate(EQ5D ~ sex + imd + age_0, data = newdata, FUN = mean)

###Code for imputing missing values for BMI with averages by age and sex

# Summarise bmi data in individuals that have data
newdata_bmi <- subset(newdata, newdata[,'bmi'] != (-1))
bmi_summary <- aggregate(bmi ~ age_0 + imd, data = newdata_bmi, FUN = mean)

# Linear model for calculating bmi - alter values of coefficients below if different
lm_bmi <- lm(bmi ~ age_0 + imd, data = bmi_summary)
summary(lm_bmi)
#plot(lm_bmi)
Intercept_bmi <- 26.761749
age_bmi <- 0.005334
imd_bmi <- 0.288154

# Calculate missing values BMI
imputed_bmi <- Intercept_bmi + newdata[,'age_0']*age_bmi + newdata[,'imd']*imd_bmi
no_bmi <- newdata[,'bmi']==(-1)
imputed_bmi <- imputed_bmi * no_bmi
newdata[,'bmi'] <- replace(newdata[,'bmi'], newdata[,'bmi']==(-1), 0)
newdata[,'bmi'] <- newdata[,'bmi'] + imputed_bmi
summary(newdata[,"bmi"])
sd(newdata[,"bmi"])

###Get table into final format and save.
##Version with all characteristics and version with only age, sex, imd and EQ5D

vars<-c("PID","sex", "age_0", "age_0","imd","ethnic","bmi","smoke","smo","units","PA","EQ5D","weighting")
#vars<-c("PID","sex","age_0","imd","EQ5D","weighting")
data<-newdata[,vars]
colnames(data)<-c("PID","sex", "age","age_0","imd","ethnic","bmi","pastsmoke","smo","units","PA","EQ5D","weighting")
#colnames(data)<-c("PID","sex","age_0","imd","EQ5D","weighting")

###Source percentiles script to get additional columns with percentiles for bmi, units and PA
#plus empty columns for family history and genetic risk.
source("HSE 2014 percentiles.R")

### Save the table
write.table(data,"X:\\ScHARR\\PR_Opt_Bowel_Cancer_Screening\\General\\PHASE 2 Optimising bowel screening\\Model\\Current CRC Model\\Data\\population.txt")
