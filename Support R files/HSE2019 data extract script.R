### CODE TO PROCESS HSE2018 DATA for EQ5D TXT FILE ###
## Lena Mandrik March 2022#

# Load up dataset

setwd("C:\\Users\\cm1om\\Documents\\Research\\Cancer\\Bladder\\model data\\HSE\\2018")

Dataset<-as.matrix(read.table("HSE2018.n.p.txt",header=TRUE))
newdata<-matrix(0,ncol=16,nrow=length(Dataset[,1]))

# Select columns to keep (column names from data dictionary)
newdata[,1]<-Dataset[,'Seriala'] # HSE serial number, Numeric
newdata[,2]<-Dataset[,'Sex'] # HSE serial number, Numeric
newdata[,3]<-Dataset[,'age16g5'] # Age categorical
newdata[,4]<-Dataset[,'qimd'] # Deprivation score, Numeric 1 = least deprived; 5 = most deprived
newdata[,5]<-Dataset[,'origin2'] # Ethnicity. Numeric 1 = white; 2 = black; 3 = Asian; 4 = mixed; 5 = other
newdata[,6]<-Dataset[,'HRPSIC7B3'] #Categorical. 3- manufacturing 
newdata[,7]<-Dataset[,'cigsta3'] # Cigarette Smoking Status. Numeric 1= current; 2 = Ex-regular smoker; 3 = Never regular smoker
newdata[,9]<-Dataset[,'Mobil17g3'] # EQ-5D response to mobility
newdata[,10]<-Dataset[,'SelfCa17g3'] # EQ-5D response to selfcare
newdata[,11]<-Dataset[,'UsualA17g3'] # EQ-5D response to Usual activity
newdata[,12]<-Dataset[,'Pain17g3'] # EQ-5D response to pain
newdata[,13]<-Dataset[,'Anxiet17g3'] # EQ-5D response to anxiety
newdata[,15]<-Dataset[,'GOR1'] # Value 3, Label = Yorkshire and the Humber
newdata[,16]<-Dataset[,'wt_int'] # Weight for analysis of core interview sample

colnames(newdata)<-c("PID","sex","age_0","imd","ethnic","occupation", "smoke", "smo", "Mobility","Selfcare","UsualAct","Pain","Anxiety","EQ5D","region", "weighting")

##Exclude data from individuals aged under 30
newdata <- subset(newdata,newdata[,'age_0']>=5)
length(newdata[,1])

#replace age to convert to continuous from categorical
age =30
category =5

for(category in 5:17){
  length.data <- nrow(subset(newdata,newdata[,'age_0']==category))
  value = round(runif(length.data, min=age, max= (age +4)), digits=0)
  newdata[ ,"age_0"] <- replace(newdata[ ,"age_0"], newdata[,'age_0']==category, values = value)
  
  category = category +1
  age = age +5
}


##Summary statistics
table(newdata[,'sex'])
summary(newdata[,'age_0'][newdata[,'age_0']>0])
sd(newdata[,'age_0'][newdata[,'age_0']>0])
table(newdata[,'age_0'])
table(newdata[,'ethnic'])
table(newdata[,'imd'])
table(newdata[,'occupation'])
table(newdata[,'region'])
table(newdata[,'smoke'])

table(newdata[,'Mobility'])
table(newdata[,'Selfcare'])
table(newdata[,'UsualAct'])
table(newdata[,'Pain'])
table(newdata[,'Anxiety'])

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

###ethnicity: assume that all those with missing data are white.
newdata[,'ethnic']<-replace(newdata[,'ethnic'],newdata[,'ethnic']==-9|newdata[,'ethnic']==-8,1)

### occupation: replace everything that is not a manufacture worker (code=3) with 0, replace 3 with 1 (binary)
newdata[,'occupation']<-replace(newdata[,'occupation'],newdata[,'occupation'] !=3, 0)
newdata[,'occupation']<-replace(newdata[,'occupation'],newdata[,'occupation'] ==3, 1)

### region: replace everything that is not Yorkshire (code=3) with 0, replace 3 with 1 (binary)
newdata[,'occupation']<-replace(newdata[,'occupation'],newdata[,'occupation'] !=3, 0)

###smoking: assume that all those with missing data were never regular smokers
newdata[,'region']<-replace(newdata[,'region'],newdata[,'region'] !=3,0)
newdata[,'region']<-replace(newdata[,'region'],newdata[,'region'] ==3,1)
table(newdata[,'region'])


###Set up Smo colomn to separate current smokers (1) from current non-smokers (0)
newdata[,'smo']<-replace(newdata[,'smo'],newdata[,'smoke']==1,1)
table(newdata[,'smo'])

###Modify smoke so that it represents past smokers
newdata[,"smoke"]<-replace(newdata[,'smoke'],newdata[,'smoke']==1|newdata[,'smoke']==3,0)
newdata[,"smoke"]<-replace(newdata[,'smoke'],newdata[,'smoke']==2,1)


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
Intercept <- 1.1042713
sex <- 0.0316358
imd <- (-0.0352753)
age <- (-0.0046997)

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


###Get table into final format and save.
##Version with all characteristics and version with only age, sex, imd and EQ5D

vars<-c("PID","sex", "age_0", "age_0","imd","ethnic","occupation","smoke","smo","region","EQ5D","weighting")


data <- newdata[,vars]
colnames(data)<-c("PID","sex", "age_0", "age","imd","ethnic","occupation","past_smoke","current_smoke","region","EQ5D","weighting")

### Save the table
write.table(data,"C:\\Users\\cm1om\\Documents\\Research\\Cancer\\Bladder\\model data\\population2018.txt")
