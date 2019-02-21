library(meta)


#######################
#this is the data
#######################
ns <- 3                                 ## number of studies
y=cbind(c(2.5,3.2,4.8),c(4,5.2,6.7))    ## mean in each arm
sd=cbind(c(0.2,0.3,0.4),c(1.4,0.9,1.1)) ## sd in eah arm
n=cbind(c(12,13,14),c(14,15,16))        ## sample size in each arm



mydata = list(ns=ns, y=y, sd=sd, n=n)


metacont(n[,1], y[,1], sd[,1], n[,2], y[,2], sd[,2],data = mydata, sm = "SMD")