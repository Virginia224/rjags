#isntallpackage
install.packages("R2jags")
install.packages("jags")
library(R2jags)
library(rjags)
help(jags)

ncont1=c(12,13,14)
ncont2=c(14,15,16)
mean1=c(187.5,193.4,161.3)
mean2=c(165.6,171.9,180.5)
sd1=c(12.3, 8.5,6.9)
sd2=c(11.9, 10.8, 11.1)
pooledSMD1=metacont(ncont1,mean1,sd1,ncont2,mean2,sd2,sm="SMD")
summary(pooledSMD1)


m=cbind(mean1,mean2)
sd=cbind(sd1,sd2)
n=cbind(ncont1,ncont2)
mydata = list(ns=3,m=m,sd=sd,n=n)

#calcualate pooled SD as data
pooled.sd.calc<-function(sd,n){
  na=length(sd)
  nominator=sqrt(apply(n * sd * sd, 1, sum))
  denominator = sqrt(apply(n, 1, sum) - na)          
pooled.sd = nominator/denominator
pooled.sd}

mydata = list(ns=3,m=m,sd=sd,n=n, pooled.sd=pooled.sd.calc(sd,n))

              
####################
#then make the model
#######################
PMAcontinuous=function() {

  for(i in 1:ns) { 
    
    #likelihood
    m[i,1] ~ dnorm(f[i,1],prec[i,1])#likelihood in one arm
    m[i,2] ~ dnorm(f[i,2],prec[i,2])#likelihood in the other arm
    
    
     f[i,1]<- u[i]*pooled.sd[i]
     f[i,2]<- (u[i] + d[i])*pooled.sd[i]
    
    #parametrisation   
    sdp[i]=((n[1,1]*sd[1,1]^2)/(n[1,1]-2))+((n[1,2]*sd[1,2]^2)/(n[1,2]-2))
    d[i] ~ dnorm(mean,prec2)
    variance[i,1]=(sd[i,1]^2)/n[i,1]
    prec[i,1]=1/variance[i,1]
    variance[i,2]=sd[i,2]^2/n[i,2]
    prec[i,2]=1/variance[i,2]
    J[i]=1-(3/(4*(n[i,1]+n[i,2])-1))
   
   }
      #prior distributions
  for (i in 1:ns) {u[i] ~ dnorm(0,.01)}
  tau ~ dunif(0,10)   #dnorm(0,100)%_%T(0,)       #This is not a proper prior for Tau!!! You are imposing FE                          
  prec2<- 1/pow(tau,2)
  mean ~ dnorm(0,0.01)

}#end of model

#####################
# initial values
#######################

initialval = NULL
#initialval = list(list(tau=0.2,mean=0.3))

#######################
# run the model
#######################

PMAinJAGS<- jags(mydata,initialval,parameters.to.save = c("tau","mean"), n.chains = 2, n.iter = 10000, n.burnin = 1000, DIC=F, model.file = PMAcontinuous)

#results
print(PMAinJAGS)