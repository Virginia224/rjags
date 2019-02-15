library(R2jags)
help(jags)


#######################
#this is the data
#######################
ns <- 3
y=cbind(c(2.5,3.2,4.8),c(4,5.2,6.7))
sd=cbind(c(0.2,0.3,0.4),c(1.4,0.9,1.1))
n=cbind(c(12,13,14),c(14,15,16))


s.within <- c() ## SD within each study
d <-  c() ## smd for each stusy
se <-  c() ##standard error of smd
for(i in 1:ns) { 
  s.within[i] <- sqrt(((n[i,1]-1)*sd[i,1]^2 + (n[i,2]-1)*sd[i,2]^2)/(n[i,1]+n[i,2]-2)) 
  d[i] <- (y[i,1] - y[i,2])/s.within[i]
  se[i] <- sqrt(((n[i,1]+n[i,2])/(n[i,1]*n[i,2])) + d[i]^2/(2*(n[i,1]+n[i,2])))
}
mydata = list(ns=ns, d=d, se=se)
#######################
#then make the model
#######################
PMAbinary=function() {
  
  for(i in 1:ns) { 
    
    #likelihood
    d[i] ~ dnorm(theta[i],prec.d[i])  #likelihood in one arm
    prec.d[i] <- 1/(se[i]^2)
        
    #parametrisation          
    theta[i] ~ dnorm(mean,prec)

    }
  
    #prior distributions
      tau ~ dunif(0,1)   #dnorm(0,100)%_%T(0,)                                 
      prec<- 1/pow(tau,2)
      mean ~ dnorm(0,100)
}#end of model

#######################
# initial values
#######################

initialval = NULL
#initialval = list(list(tau=0.2,mean=0.3))

#######################
# run the model
#######################

PMAinJAGS<- jags(mydata,initialval,parameters.to.save = c("tau","mean"), n.chains = 2, n.iter = 10000, n.burnin = 1000, DIC=F, model.file = PMAbinary)

#results
print(PMAinJAGS)

#check chain mixing
traceplot(PMAinJAGS)

#what else is there
names(PMAinJAGS)
names(PMAinJAGS$BUGSoutput)
PMAinJAGS$BUGSoutput$summary


#///////////////////////////

#######################
# parallelising
#######################

start_time <- Sys.time()
PMAinJAGS<- jags(mydata,initialval,parameters.to.save = c("tau","mean"), n.chains = 4, n.iter = 100000, n.burnin = 1000, DIC=F, model.file = PMAbinary)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
PMAinJAGS<- jags.parallel(mydata,initialval,parameters.to.save = c("tau","mean"), n.chains = 4, n.iter = 100000, n.burnin = 1000, DIC=F, model.file = PMAbinary)
end_time <- Sys.time()
end_time - start_time




