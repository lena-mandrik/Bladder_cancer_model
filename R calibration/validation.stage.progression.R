# This code uses m.BC.T.to.Stage to retrieve predicted time to other stages


# Allocate the time to stage at diagnosis for each person in HSE
m.BC.T.to.Stage <- f.stage(Mean.t.StI.StII, shape.t.StI.StII, Mean.t.StII.StIII, shape.t.StII.StIII, 
                           Mean.t.StIII.StIV, shape.t.StIII.StIV, n.i)

plot(m.BC.T.to.Stage[, "T.onsetToStage2"])
plot(m.BC.T.to.Stage[, "T.onsetToStage3"])
plot(m.BC.T.to.Stage[, "T.onsetToStage4"])

(mean.calc.T.onsetToStage2 <- mean(m.BC.T.to.Stage[ ,"T.onsetToStage2"]))
(mean.calc.T.onsetToStage3 <- mean(m.BC.T.to.Stage[ ,"T.onsetToStage3"]))
(mean.calc.T.onsetToStage4 <- mean(m.BC.T.to.Stage[ ,"T.onsetToStage4"]))

#proportion of those who progress within 1 year since onset
(pro.less1y.st2<- sum(m.BC.T.to.Stage[ ,"T.onsetToStage2"]<1)/n.i)
(pro.less1y.st3 <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage3"]<1)/n.i)
(pro.less1y.st4 <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage4"]<1)/n.i)

#proportion of those who progress within 3 years since onset
pro.less3y.st2<- sum(m.BC.T.to.Stage[ ,"T.onsetToStage2"]<3)/n.i
pro.less3y.st3 <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage3"]<3)/n.i
pro.less3y.st4 <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage4"]<3)/n.i

#proportion of those who progress more than 10 years since onset
(pro.less10y.st2<- sum(m.BC.T.to.Stage[ ,"T.onsetToStage2"]>10)/n.i)
(pro.less10y.st3 <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage3"]>10)/n.i)
(pro.less10y.st4 <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage4"]>10)/n.i)

#proportion of those who progress more than 20 years since onset
(pro.less20y.st2<- sum(m.BC.T.to.Stage[ ,"T.onsetToStage2"]>20)/n.i)
(pro.less20y.st3 <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage3"]>20)/n.i)
(pro.less20y.st4 <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage4"]>20)/n.i)

# For stage 2 proportion of patients thatprogressed between 0-7 years
(pro.st2.0.7y <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage2"]<8)/n.i)
(pro.st3.0.7y <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage3"]<12)/n.i)
(pro.st4.0.7y <- sum(m.BC.T.to.Stage[ ,"T.onsetToStage4"]<16)/n.i)

