library(R2jags)
help(jags)
library(meta)

#######################
#this is the data
#######################
ns <- 3                                 ## number of studies
y=cbind(c(2.5,3.2,4.8),c(4,5.2,6.7))    ## mean in each arm
sd=cbind(c(0.2,0.3,0.4),c(1.4,0.9,1.1)) ## sd in eah arm
n=cbind(c(12,13,14),c(14,15,16))        ## sample size in each arm
mydata = list(ns=ns, y=y, sd=sd, n=n)

# run pairwise meta-analysis using meta package
pooledSMD1=metacont(n[,1], y[,1], sd[,1], n[,2], y[,2], sd[,2],data = mydata, sm = "SMD")
summary(pooledSMD1)




#######################
#then make the model (Random effects)
#######################
PMAcontinuous=function() {
  
  for(i in 1:ns) { 
    
    #likelihood
    y[i,1] ~ dnorm(phi[i,1],prec.i[i,1])  #likelihood in one arm
    y[i,2] ~ dnorm(phi[i,2],prec.i[i,2])  #likelihood in one arm
    prec.i[i,1] <- 1/((sd[i,1]^2)/n[i,1])
    prec.i[i,2]<- 1/((sd[i,2]^2)/n[i,2])

    #calculate pooled sd
    numerator[i] <- sum(n[i,1:2]*sd[i,1:2]*sd[i,1:2])

    s.pooled[i] <- sqrt(numerator[i]/(sum(n[i,1:2])-2))

    
    #parametrisation
    phi[i,1]<- u[i]*s.pooled[i]
    phi[i,2]<- (u[i] + theta[i])*s.pooled[i]  ## SMD
    theta[i] ~ dnorm(mean,prec)
    
  }
  
  #prior distributions
  for (i in 1:ns) {
    u[i] ~ dnorm(0,.01)
  }
  tau ~ dunif(0,2)   #dnorm(0,100)%_%T(0,)                                 
  prec<- 1/pow(tau,2)
  mean ~ dnorm(0,100)
  
}
#end of model


#######################
# initial values
#######################

initialval = NULL
#initialval = list(list(tau=0.2,mean=0.3))



#######################
# run the model
#######################

PMAinJAGS<- jags(mydata,initialval,parameters.to.save = c("tau","mean", "theta"), n.chains = 2, n.iter = 10000, n.burnin = 1000, DIC=F, model.file = PMAcontinuous)

#results
print(PMAinJAGS)





#######################
# make the model (Fixed effects)
#######################
PMAcontinuous_FE=function() {
  
  for(i in 1:ns) { 
    
    #likelihood
    y[i,1] ~ dnorm(phi[i,1],prec.i[i,1])  #likelihood in one arm
    y[i,2] ~ dnorm(phi[i,2],prec.i[i,2])  #likelihood in one arm
    prec.i[i,1] <- 1/((sd[i,1]^2)/n[i,1])
    prec.i[i,2]<- 1/((sd[i,2]^2)/n[i,2])
    
    #calculate pooled sd
    numerator[i] <- sum(n[i,1:2]*sd[i,1:2]*sd[i,1:2])
    
    s.pooled[i] <- sqrt(numerator[i]/(sum(n[i,1:2])-2))
    
    
    #parametrisation
    phi[i,1]<- u[i]*s.pooled[i]
    phi[i,2]<- (u[i] + mean)*s.pooled[i]  ## SMD
    # theta[i] ~ dnorm(mean,prec)
    
  }
  
  #prior distributions
  for (i in 1:ns) {
    u[i] ~ dnorm(0,.01)
  }
  tau ~ dunif(0,2)   #dnorm(0,100)%_%T(0,)                                 
  prec<- 1/pow(tau,2)
  mean ~ dnorm(0,100)
  
}
#end of model


#######################
# initial values
#######################

initialval = NULL
#initialval = list(list(tau=0.2,mean=0.3))



#######################
# run the model
#######################

PMAinJAGS_FE<- jags(mydata,initialval,parameters.to.save = c("mean"), n.chains = 2, n.iter = 10000, n.burnin = 1000, DIC=F, model.file = PMAcontinuous_FE)

#results
print(PMAinJAGS_FE)