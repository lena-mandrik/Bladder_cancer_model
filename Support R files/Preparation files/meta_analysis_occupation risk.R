# metafor meta analysis RR occupation

library(metafor)

dat.bcg

dat <- escalc(measure="RR", ai = tpos, bi=tneg,
              ci=cpos, di=cneg,
              slab = paste(author, ",", year, sep=""), data=dat.bcg)


m.dat <- cbind(c(1,2,3), c(1969,1970,1976), c(10,1,37), c(4,1,19), c(1400,2120, 40876),c(1400,2120, 40876))
colnames(m.dat) <- c("author", "year", "cases", "expected", "N_cases", "N_control")
m.dat <- (as.data.frame(m.dat))

as.numeric(m.dat[ ,3:6])

m.dat[ ,1] <- c("Veys","Anthony","Fox")

dat <- escalc(measure="RR", ai = cases, bi=N_cases,
              ci=expected, di=N_control,
              slab = paste(author, ",", year, sep=""), data=m.dat)

dat

# random-effect model using log risk ratios and variances as input
res <- rma(yi, vi, data=dat)
res

# Predicted pooled risk ratio and corresponding CI
predict(res, transf=exp, digits =2)

#forest plot
forest(res)

forest(res, addpred =TRUE, header =T)
