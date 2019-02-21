
#*********************************************************************************
#             Load the libraries needed             
#*********************************************************************************
library(meta)
library(metafor)
library(netmeta)
library(readxl)

library(devtools)
library(remotes)
install_github("esm-ispm-unibe-ch/NMAJags")
library(NMAJags)
library(R2jags)


#*********************************************************************************
#             NMA USING JAGS              
#*********************************************************************************


Griselda <- read_excel("Depression.xlsx")

ls("package:NMAJags")

make.jagsNMA.data

table(Griselda$drug_name)
#transform the data into a list suitable for JAGS analysis
NMAdataBinary=make.jagsNMA.data(studyid=studyID,t=drug_name,r=Responders,n=Ntotal,data=Griselda,type="binary",reference = "Amitriptyline")

modelNMABinary

#run Jags and create a jags object
NMAinJAGS<- jags(data = NMAdataBinary, inits = NULL,
                 parameters.to.save = c("ORref","tau"), n.chains = 2, n.iter = 10000,
                 n.burnin = 1000,DIC=F,n.thin=10,
                 model.file = modelNMABinary)
print(NMAinJAGS)
save(NMAinJAGS,file="NMAinJAGSall.RData",envir = .GlobalEnv)

#check chain mixing
traceplot(NMAinJAGS,varname="tau" )
traceplot(NMAinJAGS,varname="ORref" )

#forestplot against placebo
y=NMAinJAGS$BUGSoutput$mean$ORref
ES=y
seES=NMAinJAGS$BUGSoutput$sd$ORref
studlab <- factor(sort(unique(Griselda$drug_name)))
studlab <- studlab[-2]
Griselda_forest <- metagen(ES, seES, studlab = studlab)
forest(Griselda_forest)



#*********************************************************************************
#             NMR USING JAGS              
#*********************************************************************************

#create random variable "year of randomisation"
rand.year <- c()
rand.year[1] <- sample(1998:2010, 1)
for (i in 2:nrow(Griselda)) {
  if(Griselda$studyID[i]==Griselda$studyID[i-1]) {
    rand.year[i] <- rand.year[i-1]
  }
  else {
    rand.year[i]=sample(1998:2010, 1, replace = T)
  }
}


#transform the data into a list suitable for JAGS analysis; this including the "year of randomisaion" variable
NMRdataBinary=make.jagsNMA.data(studyid=studyID,t=drug_name,r=Responders,n=Ntotal,data=Griselda,type="binary",reference = "Amitriptyline", othervar = rand.year)

# define NMR model for binary outcome
modelNMRbinary <- function () 
{
  for (i in 1:ns) {
    w[i, 1] <- 0
    theta[i, t[i, 1]] <- 0
    for (k in 1:na[i]) {
      r[i, t[i, k]] ~ dbin(p[i, t[i, k]], n[i, t[i, k]])
    }
    logit(p[i, t[i, 1]]) <- u[i]
    for (k in 2:na[i]) {
      logit(p[i, t[i, k]]) <- u[i] + theta1[i, t[i, k]]
      theta1[i, t[i, k]] <- theta[i, t[i, k]] + beta[t[i, 1], t[i, k]] * variab[i]
      theta[i, t[i, k]] ~ dnorm(md[i, t[i, k]], precd[i, t[i, k]])
      md[i, t[i, k]] <- mean[i, k] + sw[i, k]
      w[i, k] <- (theta[i, t[i, k]] - mean[i, k])
      sw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)
      precd[i, t[i, k]] <- prec * 2 * (k - 1)/k
      mean[i, k] <- d[t[i, k]] - d[t[i, 1]]
    }
  }
  for (i in 1:ns) {
    u[i] ~ dnorm(0, 0.01)
  }
  tau ~ dnorm(0, 1) %_% T(0, )
  prec <- 1/pow(tau, 2)
  tau.sq <- pow(tau, 2)
  d[ref] <- 0
  for (k in 1:(ref - 1)) {
    d[k] ~ dnorm(0, 0.01)
  }
  for (k in (ref + 1):nt) {
    d[k] ~ dnorm(0, 0.01)
  }
  for (i in 1:(nt - 1)) {
    for (j in (i + 1):nt) {
      OR[j, i] <- exp(d[j] - d[i])
      LOR[j, i] <- d[j] - d[i]
    }
  }
  for (j in 1:(ref - 1)) {
    ORref[j] <- exp(d[j] - d[ref])
  }
  for (j in (ref + 1):nt) {
    ORref[j] <- exp(d[j] - d[ref])
  }
  for (i in 1:nt) {
    for (j in 1:nt) {
      beta[i, j] <- b[j] - b[i]
    }
  }
  b[ref] <- 0
  for (k in 1:(ref - 1)) {
    b[k] ~ dnorm(0, 1e-04)
  }
  for (k in (ref + 1):nt) {
    b[k] ~ dnorm(0, 1e-04)
  }
}

#run Jags and create a jags object
NMRinJAGS<- jags(data = NMRdataBinary, inits = NULL,
                 parameters.to.save = c("ORref","tau", "beta"), n.chains = 2, n.iter = 10000,
                 n.burnin = 1000,DIC=F,n.thin=10,
                 model.file = modelNMRbinary)
print(NMRinJAGS)
save(NMRinJAGS,file="NMRinJAGSall.RData",envir = .GlobalEnv)
