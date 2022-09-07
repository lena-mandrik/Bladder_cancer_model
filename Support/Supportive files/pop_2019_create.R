setwd("C:\\Users\\cm1om\\Documents\\Research\\Cancer\\Bladder")

pop19 <- read.table("model data\\HSE19.txt", header = T,sep="\t")

pop14 <- read.table("model data\\HSE14.txt", header = T)


EQ5D <- rep(0, nrow(pop19))

pop19 <-cbind(pop19, EQ5D)

i=0
n=1
m=1
v=30

for(i in 0:1){ #sex
  for(n in 1:5){ #imd
    for(m in 1:5){ #ethnic
     for(v in 30:100){ #age
       
       condition <- sum(pop14[ , "sex"]== i & pop14[ , "imd"]== n & pop14[ , "ethnic"]== m & pop14[ , "age"] >= (v-5) & pop14[ , "age"] <= (v+5))
      
       #for(l in 1:nrow(pop19)){
        size_pop <-sum(pop19[, "sex"]== i & pop19[, "imd"]== n & pop19[, "ethnic"]== m & pop19[, "age"] >= (v-5) & pop19[, "age"] <= (v+5))
         
         if(condition>0 & size_pop >0){
           pop19[pop19[, "sex"]== i & pop19[, "imd"]== n & pop19[, "ethnic"]== m & pop19[, "age"] >= (v-5) & pop19[, "age"] <= (v+5),"EQ5D"] <- sample(pop14[pop14[ , "sex"]== i & pop14[ , "imd"]== n & pop14[ , "ethnic"]== m & pop14[ , "age"] >= (v-5) & pop14[ , "age"] <= (v+5), "EQ5D"],size_pop,replace=T)
           
         }       
           v=v+5
           
       }
        m=m+1
        
     }
      n=n+1
      
    } 
    i = i+1
}

mean(pop19[ ,"EQ5D"])
mean(pop14[ ,"EQ5D"])

pop19[ ,"EQ5D"] <- replace(pop19[ ,"EQ5D"], pop19[ ,"EQ5D"]==0, median(pop14[ ,"EQ5D"]))
  
write.table(pop19, file= "model data\\population.txt", col.names= T, sep="\t") 

read.table("model data\\population.txt", header = T,sep="\t")
