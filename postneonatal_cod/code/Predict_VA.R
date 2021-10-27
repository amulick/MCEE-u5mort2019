
########################################
###
###   Analysis of final data including two periods for separate indian states
###
###   David Prieto (Feb 2020)
###
### Top ten burden countries (1-59 months)
###
### Nigeria (NGA)
### India (IND) ---- envelopes, aids and measles?
### DRC (COD)
### Pakistan (PAK)
### Ethiopia (ETH)
### China (CHN)
### Tanzania (TZA)
### Indonesian (IDN)
### Angola (AGO)
### Bangladesh (BGD)
### Niger (NER)
###
###
########################################

#####################################################################################################################

## China
## India

## PostNeonates
##   VR modeled
##   VA modeled
##      -Bayes
##      -post hoc
##   Good VR

## Measles
## HIV
## WHO malaria
## Envelopes

#####################################################################################################################

rm(list=ls())
Start<- Sys.time()



#### Analysis of samples
library(reshape2)
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
#load(file="../code_Jags_2020_08 31 wo INCAM/PN_data_2020_08_31_INCAM.RData")
load(file="../data/PN_VA_data.RData")
vxf
summary(studies)
names(studies)

#for NIUE 2005
data.predict$envelope_aids_measles_free<- ifelse(data.predict$envelope_aids_measles_free==0,1,
                                                 data.predict$envelope_aids_measles_free)
dim(data.predict)
table(data.predict$year)
data.predict<- data.predict[which(data.predict$year>1999),]
table(data.predict$year)
length(unique(data.predict$isocode))
summary(data.predict$envelope_aids_measles_free==0)


#for India States
#table(data.predict$isocode)


###################################X
#######   PREDICTIONS     ##########
###################################X
source("functions_predictVA.R")  ##  functions

#Contains object named "out"
#load(file="../code_Jags_2020_08 31 wo INCAM/results_VA/f.e1_Test24d_INCAM_lambda_10_sd_0.14.RData")
#out$param$VXF
#out$output$BUGSoutput$mean$B
load(file="../results/VA.RData")
out$output$BUGSoutput$mean$B
vxf<- out$param$VXF
vxf




#check ranges/scales
for(i in out$param$VXF){
  print(i)
  print(summary(studies[,i]))
  print(summary(data.predict[,i]))
}

#convergence summary
summ<- as.data.frame(out$output$BUGSoutput$summary)
summary(summ$Rhat)


#function(MO, STUP, NP=500, IDV=c("isocode","year"), PPR=c(0, 0))
#iterations saved ,  chains, parameters
dim(out$output$BUGSoutput$sims.array)
p21 <- f.ci2(out, data.predict, nset=200, id=c("isocode","year"),PRAN=T)
dim(p21[[1]])
dim(data.predict)
length(unique(data.predict$isocode))
#studies[which(studies[,vxr] %in% data.predict[,vxr]),vxr]
#studies[which(studies[,vxr] %in% data.predict[,vxr]),"study_id"]
dim(p21[[2]])
#rownames(p21[[1]])

####################################
#######   WHO REGIONS   ############
####################################
codes<- read.dta13("../../other_inputs/codelist region.dta")
#head(codes)
codes$isocode<- codes$iso3


####################################
#######   WHO MALARIA   ############
####################################
data.predict$iso.nat<- substr(data.predict$isocode,1,3)
#table(data.predict$iso.nat)
dp.nat<- aggregate(data.predict[,c("envelope_aids_measles_free","pfpr")],
                   by=list(data.predict$year,
                           data.predict$iso.nat), FUN=sum)
names(dp.nat)<- c("year","isocode","envelope_aids_measles_free","pfpr")
#head(dp.nat)
dp.nat$pfpr<- ifelse(dp.nat$isocode=="IND", dp.nat$pfpr/21, dp.nat$pfpr)
#head(dp.nat)

malunc<- read.dta("../data/malaria_2000_2019.dta")
malunc<- malunc[which(malunc$iso3 %in% dp.nat$isocode),]


rownames(malunc)<- paste(malunc$iso3,malunc$year,sep=".")
library(arrayhelpers)
source("posthoc_hib.R")
source("addWHOmal.R")
#names(data.predict)
#head(data.predict[,c("isocode","year","envelope")])










NSim<- dim(p21[[2]])[3]

library(plyr)
#data.predict$hib3

#names(data.predict)
############################################################
#Post hoc adjustment

dim(p21$estR)
head(p21$estR[,,1])


xt<- which(data.predict$pfpr==0)
summary(p21$estR[xt,,1])

#for fixed and random
test3<-apply(p21$estR, 3, FUN=posthoc_hib, verbose=FALSE)
#for fixed only
#test3<-apply(p21$estF, 3, FUN=posthoc)

#Transform list back to array
testPH<-array(unlist(test3), dim = c(nrow(test3[[1]]), ncol(test3[[1]]), length(test3)))
dimnames(testPH)[[1]]<-rownames(p21$estR)
dimnames(testPH)[[2]]<-colnames(p21$estR)
p21.PH<- p21
p21.PH$estR<- testPH
dim(p21.PH$estR)
#summary(p21.PH$estR[,,1])

############################################################
#Aggregate to national India
source("aggregate_India.R")
test4<-apply(p21.PH$estR, 3, FUN=aggregate_India)
#Transform list back to array
testIND<-array(unlist(test4), dim = c(nrow(test4[[1]]), ncol(test4[[1]]), length(test4)))


summary(data.predict$envelope_aids_measles_free)

dimnames(testIND)[[1]]<-  paste0(dp.nat$isocode,".", dp.nat$year)
dimnames(testIND)[[2]]<-colnames(p21$estR)
p21.IND<- p21
p21.IND$estR<- testIND
dim(p21.IND$estR)
head(p21.IND$estR[,,1])


############################################################
#WHO malaria
source("addWHOmal.R")
test2<-apply(p21.IND$estR, 3, FUN=addWHOmal)
warnings()
vxf
#Transform list back to array
testA<-array(unlist(test2), dim = c(nrow(test2[[1]]), ncol(test2[[1]]), length(test2)))
dim(testA)
dimnames(testA)[[1]]<-paste0(dp.nat$isocode,".", dp.nat$year)
dimnames(testA)[[2]]<-colnames(p21$estR)
p21.whoMAL<- p21
p21.whoMAL$estR<- testA


#Make data set of country point estimates
dim(testPH)
est<-apply(testA, c(1,2), FUN=mean)
lower<-apply(testA, c(1,2), FUN=quantile, probs=c(0.025))
upper<-apply(testA, c(1,2), FUN=quantile, probs=c(0.975))
colnames(est) <-colnames(lower)<-colnames(upper)<-vdr
rownames(est) <-rownames(lower)<-rownames(upper)<-paste0(dp.nat$isocode,".", dp.nat$year) #rownames(p21$estR)
colnames(lower)<- paste(colnames(lower),".lower",sep="")
colnames(upper)<- paste(colnames(upper),".upper",sep="")


VA.est<- merge(est, lower, by="row.names")
rownames(VA.est)<- VA.est$Row.names
VA.est$Row.names<- NULL
VA.est<- merge(VA.est, upper, by="row.names")
#rownames(VA.est)<- VA.est$Row.names
temp<-strsplit(rownames(VA.est), split="\\.")
head(temp)
VA.est$iso3<- dp.nat$isocode
VA.est$year <- dp.nat$year
head(VA.est)
dim(VA.est)
VA.est<- merge(VA.est, codes[,c("iso3","whoreg6")], by="iso3")
head(VA.est)
dim(VA.est)

VA.est$Row.names<- NULL
head(VA.est)
VA.draws<- testA
VA.draws.noPH<- p21$estR












# To aggregate
for( i in vdt) {
  VA.est[,i]<- VA.est[,i]* dp.nat$envelope_aids_measles_free
}

table(VA.est$year)
dim(VA.est)
agg.new<-aggregate(VA.est[, c(vdt)], list(VA.est$year), FUN=sum)
agg.new$tot<- rowSums(agg.new[,vdt])
length(unique(dp.nat$isocode))
length(unique(dp.nat$year))



VA.est.noPH<-as.data.frame(apply(p21$estR, c(1,2), FUN=mean))
head(VA.est.noPH)
rownames(VA.est.noPH)<- rownames(p21$estR)
colnames(VA.est.noPH)<- vdr
temp<-strsplit(rownames(VA.est.noPH), split="\\.")
#EVA.est.noPH$year<-unlist(lapply(temp, FUN=function(x) x[2]))
#VA.est.noPH$iso3<-unlist(lapply(temp, FUN=function(x) x[1]))
VA.est.noPH$iso3<- data.predict$isocode
VA.est.noPH$year <- data.predict$year

dim(VA.est)
dim(VA.est.noPH)
#agg.newM
for( i in vdt) {
  VA.est.noPH[,i]<- VA.est.noPH[,i]* data.predict$envelope_aids_measles_free
}
VA.est.noPH$year<- data.predict$year
VA.est.noPH$iso3<- data.predict$isocode

VA.est.noPH<- merge(VA.est.noPH, codes[,c("iso3","whoreg6")], by="iso3")
head(VA.est.noPH)
dim(VA.est.noPH)

save(VA.draws, VA.est,
     file=paste0("../results/Predict_VA.RData"))





#COEFFICIENTS
coef<- data.frame( cause = rep(vdt[-1],each=np), cause.n=rep(2:8, each=np),
                   cov=rep(c("Intercept",vxf), 7), est=as.vector(out$output$BUGSoutput$mean$B))
head(coef)
############################################


#VA.est$pnd<- rowSums(VA.est[,vdt])
q("no")















