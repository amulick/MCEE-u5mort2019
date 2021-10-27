
########################################Z
###
###   Analysis of final data including two periods for separate indian states
###
###   David Prieto (Feb 2020)
###
########################################

rm(list=ls())
library(coda)
library(lattice)
library(foreign)
library(tidyverse)
library(ggpubr)
library(abind)
library(readstata13)



#Predictions from VA model
load("../results/Predict_VA.RData")



#Predictions from VR model
load("../results/Predict_VR.RData")
dimnames(VR.draws)[[2]]
length(unique(VR.est$iso3))
names(VR.est)<- gsub("perinatal","neonatal", names(VR.est))
dimnames(VR.draws)[[2]]<- gsub("perinatal","neonatal", dimnames(VR.draws)[[2]])

unique(VA.est$iso3) [!unique(VA.est$iso3) %in% unique(VR.est$iso3)]


# load a VA and VR models to get variables needed for prediction
# VR model
load("../results/VR.RData")
out$param$VXF
out$param$VDT
vr.out<- out
rownames(VA.draws)
it<- rownames(VA.draws)
it<- substr(it, 1,3)

#VA model
load("../results/VA.RData")
va.out<- out
va.out$param$VXF
va.out$param$VDT


## Load data from all countries:
dc <- read.dta13(file="../data/master national level data_1 May 2020_1980-2019.dta")
dc$u5mr<- dc[,"_5q0"]
summary(dc$u5mr)
table(dc$iso3)
dp<- dc[which(dc$iso3 %in% it & dc$year>=2000),]
## di <- read_csv(file="india_covs_1990-2019_1May2020.csv")

###  Combine VA and VR predictions
###
## For each country/year crate a weight of VA dep
table(dp$iso3)
dp$isoyear <- paste(dp$iso3, dp$year, sep=".")
dp$wei.vavr <- sapply(dp$u5mr, function(x) max(0, min(1, (x - 25)/ (35 - 25) )))


# causes in each model
ca <- va.out$param$VDT

cr <- vr.out$param$VDT
cr<- gsub("perinatal","neonatal", cr)
cc <- unique(c(cr,ca))




## Complete matrices with non-existing causes
# in VA matrices
dim(VA.draws)
PFa <- abind(VA.draws, array(0, dim=c(dim(VA.draws)[1],
                                      length(cc[!(cc %in% ca)]), dim(VA.draws)[3])), along=2)
dimnames(PFa)[[2]] <- c(ca, cc[!(cc %in% ca)])
dimnames(PFa)[[1]]
PFa <- PFa[dp$isoyear,cc,]
PFa[1:3,,1:2]


PRa <- abind(VA.draws, array(0, dim=c(dim(VA.draws)[1],
                                      length(cc[!(cc %in% ca)]), dim(VA.draws)[3])), along=2)
dimnames(PRa)[[2]] <- c(ca, cc[!(cc %in% ca)])
PRa <- PRa[dp$isoyear,cc,]
PRa[1:3,,1:2]


# in VR matrices
PFr <- abind(VR.draws, array(0, dim=c(dim(VR.draws)[1],
                                      length(cc[!(cc %in% cr)]), dim(VR.draws)[3])), along=2)
dimnames(PFr)[[2]] <- c(cr, cc[!(cc %in% cr)])
PFr <- PFr[dp$isoyear,cc,]

table( dp$isoyear %in% dimnames(PFr)[[1]] )
dp$isoyear[which ( !dp$isoyear %in% dimnames(PFr)[[1]] )]



PFr[1:3,,1:2]
PRr <- abind(VR.draws, array(0, dim=c(dim(VR.draws)[1],
                                      length(cc[!(cc %in% cr)]), dim(VR.draws)[3])), along=2)
dimnames(PRr)[[2]] <- c(cr, cc[!(cc %in% cr)])
PRr <- PRr[dp$isoyear,cc,]


## array of weights
WM <- array(dp$wei.vavr, dim=dim(PFa))

xt<-dimnames(PFa[which(!WM[,1,5] %in% c(0,1)),,1])[[1]]

## Combination
pco.a <- list(Model="VAVR",
              PF=(PFa*WM + PFr*(1-WM)), PR=(PRa*WM + PRr*(1-WM)))
MA.draws<- pco.a
# Calculate predictions
est<-apply(pco.a$PF, c(1,2), FUN=mean)



save(est, MA.draws, file="../results/Model_averaging.RData")
#write.csv(est, file="Prediction2_model_averaging_JP.csv")
#write.dta(as.data.frame(est), file="../results/Prediction2_model_averaging_JP.dta")

q('no')
