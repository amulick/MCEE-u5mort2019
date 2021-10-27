
########################################
###
###   Analysis of final data including two periods for separate indian states
###
###   David Prieto (Feb 2020)
###
###
###
########################################


rm(list=ls())
Start<- Sys.time()



#### Analysis of samples
library(R2jags)
library(doParallel)
library(coda)
library(lattice)
library(foreign)
#library(saveJAGS)
library(readstata13)

###
###  Prepare datasets for R
###

# Read data
load(file="../data/PN_VR_data.RData")
dim(data.predict)
vxf #covariates
vdt #causes
dim(studies)
table(studies$year)

###################################X
#######   PREDICTIONS     ##########
###################################X
source("functions_predictVR.R")  ##  functions

#Contains object named "out"
#load(file="../results/VR_Bayes_lambda5_sd_0.07.RData")
load(file="../results/VR.RData")

ls()
summary(out$Studies$year)
dim(out$output$BUGSoutput$sims.array)
dimnames(out$output$BUGSoutput$sims.array)

vxf<- out$param$VXF
vxf


dim(out$output$BUGSoutput$sims.array)
dim(data.predict)




#data.predict<- data.predict[which(data.predict$isocode %in% prev$iso3),]
data.predict<- data.predict[which(#data.predict$isocode !="IND" &
                                  data.predict$year !=2020
                                  & data.predict$year>1999
                                  ),]


#iterations saved, chains, parameters
dim(out$output$BUGSoutput$sims.array)
p21 <- f.ci2(out, data.predict, nset=200, id=c("isocode","year"), PRAN=T)


####################################
#######   WHO REGIONS   ############
####################################
codes<- read.dta13("../../other_inputs/codelist region.dta")
head(codes)
codes$isocode<- codes$iso3


library(arrayhelpers)


NSim<- dim(p21[[2]])[3]
library(plyr)


#Make data set of country point estimates
dim(p21$estR)
est<-apply(p21$estR, c(1,2), FUN=mean, na.rm=TRUE)
lower<-apply(p21$estR, c(1,2), FUN=quantile, probs=c(0.025), na.rm=TRUE)
upper<-apply(p21$estR, c(1,2), FUN=quantile, probs=c(0.975), na.rm=TRUE)
colnames(est) <-colnames(lower)<-colnames(upper)<-vdr
rownames(est) <-rownames(lower)<-rownames(upper)<-rownames(p21$estR)
colnames(lower)<- paste(colnames(lower),".lower",sep="")
colnames(upper)<- paste(colnames(upper),".upper",sep="")

dim(est)
test<- rowSums(est)
summary(test)
head(test)

VR.est<- merge(est, lower, by="row.names")
rownames(VR.est)<- VR.est$Row.names
VR.est$Row.names<- NULL
VR.est<- merge(VR.est, upper, by="row.names")
#rownames(VR.est)<- VR.est$Row.names
temp<-strsplit(rownames(VR.est), split="\\.")

VR.est$iso3<- data.predict$isocode
VR.est$year <- data.predict$year
VR.est<- merge(VR.est, codes[,c("iso3","whoreg6")], by="iso3")

VR.est$Row.names<- NULL
head(VR.est)
VR.draws<- p21$estR


#dounble check intervals
for(C in vdr) {
  print(summary(VR.est[,paste0(C,"")] - VR.est[,paste0(C,".lower")]))
}

summary(VR.est)
summary(VR.est[,vdt])

save(VR.draws, VR.est, file=paste0("../results/Predict_VR.RData"))



#COEFFICIENTS
np<- length(vxf)+1
coef<- data.frame( cause = rep(vdt[-1],each=np), cause.n=rep(2:8, each=np),
                   cov=rep(c("Intercept",vxf), 7), est=as.vector(out$output$BUGSoutput$mean$B))
head(coef)
q("no")
