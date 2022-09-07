# this script calculates characteristic for different populations 


setwd("C:\\Users\\cm1om\\Documents\\Research\\Cancer\\Bladder")

pop <- read.table("model data\\population.txt", header = T,sep="\t")


# Create representative population to England

Probability <- pop[ ,"weight"]/10

pop <-cbind(pop, Probability)

pop_new <- pop[1,]

i =1

for(i in 1:100){
  
  n_row_pop <- sample(1:nrow(pop), size=1000, replace=F, prob = pop[ ,"Probability"])
  
  pop_test <- pop[n_row_pop,]
  
  pop_new <- rbind(pop_new, pop_test)
  
  i = i +1
}


mean(pop[ ,"age"])
sd(pop[ ,"age"])
mean(pop_new[ ,"age"])
sd(pop_new[ ,"age"])

mean(pop[ ,"sex"])
sum(pop[ ,"sex"])
mean(pop_new[ ,"sex"])


mean(pop[ ,"sex"])


sum(pop_new[ ,"ethnic"]==1)/nrow(pop_new)
sum(pop_new[ ,"ethnic"]==2)/nrow(pop_new)
sum(pop_new[ ,"ethnic"]==3)/nrow(pop_new)

sum(pop_new[ ,"imd"]==1)/nrow(pop_new)
sum(pop_new[ ,"imd"]==2)/nrow(pop_new)
sum(pop_new[ ,"imd"]==3)/nrow(pop_new)
sum(pop_new[ ,"imd"]==4)/nrow(pop_new)
sum(pop_new[ ,"imd"]==5)/nrow(pop_new)

sum(pop_new[ ,"smoker"])/nrow(pop_new)
sum(pop_new[ ,"past_smoker"])/nrow(pop_new)


sum(pop[ ,"past_smoker"]==1)

mean(pop_new[ ,"sex"])

# Sample the population for Yorkshire
pop_Y <- pop_new[pop_new[ ,"region"]==1, ]


sum(pop[ ,"occupation"])


sum(pop_Y[ ,"ethnic"]==1)/nrow(pop_Y)
sum(pop_Y[ ,"ethnic"]==2)/nrow(pop_Y)
sum(pop_Y[ ,"ethnic"]==3)/nrow(pop_Y)

sum(pop_Y[ ,"imd"]==1)/nrow(pop_Y)
sum(pop_Y[ ,"imd"]==2)/nrow(pop_Y)
sum(pop_Y[ ,"imd"]==3)/nrow(pop_Y)
sum(pop_Y[ ,"imd"]==4)/nrow(pop_Y)
sum(pop_Y[ ,"imd"]==5)/nrow(pop_Y)

sum(pop_Y[ ,"smoker"])/nrow(pop_Y)
sum(pop_Y[ ,"past_smoker"])/nrow(pop_Y)
