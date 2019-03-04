
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
forest(Griselda_forest, Griselda_forest$TE, overall=F) ###sorted by Treatment effect



#*********************************************************************************
#             NMR USING JAGS              
#*********************************************************************************

#create random variable "year of randomisation" ###!!!!!??? why????/
#rand.year <- c()
#rand.year[1] <- sample(1998:2010, 1)
#for (i in 2:nrow(Griselda)) {
#  if(Griselda$studyID[i]==Griselda$studyID[i-1]) {
#    rand.year[i] <- rand.year[i-1]
#  }
#  else {
#    rand.year[i]=sample(1998:2010, 1, replace = T)
#  }
#}
#Griselda$rand.year.c=rand.year-2004

### simpler code to generate random variable "year of randomisation"
a=table(Griselda$studyID)
rand.year.c=sample(-5:5,length(a),replace=T)
Griselda$rand.year.c=rep(rand.year.c,a)
####

#transform the data into a list suitable for JAGS analysis; including the "year of randomisation" variable
NMRdataBinary=make.jagsNMA.data(studyid=studyID,t=drug_name,r=Responders,n=Ntotal,data=Griselda,type="binary",reference = "Amitriptyline", othervar = rand.year.c)

source("modelNMRbinary.R")


#run Jags and create a jags object
NMRinJAGS<- jags(data = NMRdataBinary, inits = NULL,
                 parameters.to.save = c("ORref","tau", "b"), n.chains = 2, n.iter = 20000,
                 n.burnin = 1000,DIC=F,n.thin=5,
                 model.file = modelNMRbinary)
print(NMRinJAGS)
save(NMRinJAGS,file="NMRinJAGSall.RData",envir = .GlobalEnv)


# check convergence####
traceplot(NMRinJAGS, varname="tau")
traceplot(NMRinJAGS, varname="ORref")
traceplot(NMRinJAGS, varname="b")

###!!!!! COMPARE THE TAU IN META-ANALYSIS WITH TAU IN REGRESSION TO SEE DIFFERENCES
print(NMAinJAGS$BUGSoutput$mean$tau)
print(NMRinJAGS$BUGSoutput$mean$tau)


