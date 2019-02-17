library(R2jags)
help(jags)


#######################
#this is the data
#######################
ns <- 3                                 ## number of studies
y=cbind(c(2.5,3.2,4.8),c(4,5.2,6.7))    ## mean in each arm
sd=cbind(c(0.2,0.3,0.4),c(1.4,0.9,1.1)) ## sd in eah arm
n=cbind(c(12,13,14),c(14,15,16))        ## sample size in each arm



mydata = list(ns=ns, y=y, sd=sd, n=n)


#######################
#then make the model (Random effects)
#######################
PMAcontinuous=function() {
  
  for(i in 1:ns) { 
    
    #likelihood
    y[i,1] ~ dnorm(phi[i,1],var.d[i,1])  #likelihood in one arm
    y[i,2] ~ dnorm(phi[i,2],var.d[i,2])  #likelihood in one arm
    var.d[i,1] <- (sd[i,1]^2)/n[i,1]
    var.d[i,2]<- (sd[i,2]^2)/n[i,2]

    s.pooled[i] <- sqrt((n[i,1]*sd[i,1]^2)/(n[i,1]-2) + (n[i,2]*sd[i,2]^2)/(n[i,2]-2))
    J[i] <- 1 - 3/(4*(n[i,1]+n[i,2])-1)
    
    
    #parametrisation
    phi[i,1] <- (u[i]*s.pooled[i])/J[i]
    phi[i,2]<- (u[i] + theta[i])*s.pooled[i]/J[i]  ## SMD
    theta[i] ~ dnorm(mean,prec)
    
  }
  
  #prior distributions
  for (i in 1:ns) {
    u[i] ~ dnorm(0,.01)
  }
  tau ~ dunif(0,1)   #dnorm(0,100)%_%T(0,)                                 
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