
########################################Z
###
###   Analysis of final data including two periods for separate indian states
###
###   David Prieto (Feb 2020)
###
########################################Z


#### Analysis of samples
rm(list=ls())
library(coda)
library(lattice)
library(foreign)
library(tidyverse)
library(ggpubr)
library(abind)
library(readstata13)


source("Prediction_functions.R")  ##  functions


label<- "20200915"

## Load data from all countries:
#dc<- read.dta13("/Users/jamieperin/Documents/Projects/COD/COD_2020/Other_inputs/master national level data_1 Sep 2020_1980-2019.dta")
dc<- read.dta13("../data/master national level data_1 Sep 2020_1980-2019.dta")
table(dc$iso3, useNA="always")

# Load model data
load(file="../data/data_missing_2020_05_27.RData")



# load a VA and VR models to get variables needed for prediction
load("../results/VA_225_0.14.RData")
vxf<-va_225_0.14$param$VXF
vxf

for(v in vxf[3:length(vxf)]){
  cat("\n\n\n",v,"\n")
print(summary(studies[,v]))
print(summary(data.predict[,v]))
}

va_225_0.14$param$VDT #causes
dim(va_225_0.14$Studies)
names(va_225_0.14$Studies)
length(unique(va_225_0.14$Studies$reterm))

load("../results/VR1.RData")
vr1_325$param$VXF
vr1_325$param$VDT #causes
head(vr1_325$Studies)
length(unique(vr1_325$Studies$isocode))
load("../results/VR2.RData")
vr2_400$param$VXF
vr2_400$param$VDT #causes
dim(vr2_400$Studies)
length(vr2_400$Deaths)
names(vr1_325)
names(vr1_325$Studies)
vr1_325$output$n.iter




# Covariates for prediction for both models
vpr <- unique(c(va_225_0.14$Vars$xvar, vr1_325$Vars$xvar, vr2_400$Vars$xvar))


# Prepare dataset for prediction
summary(dc$iso3)
names(dc) <- gsub("_5q0", "u5mr",names(dc))
head(dc[which(is.na(dc$sba)),c("iso3","year")])


dp <- transmute(dc, isocode=iso3, country=whoname, year=year, u5mr=u5mr_sm, nmr=nmr_sm, dpt=dtp3_sm, gfr=gfr_sm,
                lbwrate=lbw_sm, sba=sba_sm, bcg=bcg_sm, pab=pab_sm, femlit=literacy_fem_sm, gni=gni_sm, gini=gini_sm, anc=anc4_sm, premvslbw=1)

## Remove countries with high quality data
dvh <- read_csv(file="../data/VR_highqual_mostrecent.csv")
dp <- filter(dp, !(isocode %in% dvh$isocode))

## put dummies for for SSA and SA
lc <- readxl::read_xlsx("../data/WPP2019_F01_LOCATIONS.XLSX", sheet="DP")
lc$regSSA <- as.numeric(lc$region == "Sub-Saharan Africa")
lc$regSA  <- as.numeric(lc$subregion == "Southern Asia")
dp <- merge(dp, lc[,c("isocode","regSA", "regSSA")], all.x=T, all.y=F)


table(dp$year)


### Add data from India (use _vsm covariates)
di <- read_csv(file="../data/india_covs_1990-2019_1May2020.csv") # covariates
di <- transmute(di, isocode=paste0("I",substr(di$state,1,3)), country=state, year=year, u5mr=u5mr, nmr=nmr, dpt=dpt_vsm, lbwrate=lbw_vsm, sba=sba_vsm, bcg=bcg_vsm, anc=anc4_vsm, femlit=lit_vsm, premvslbw=1)
## The Indian states don´t have: gfr, pab=pab_sm, gni=gni_sm, gini=gini_sm,
# Add it from the data at national level
di <- left_join(di, dp[dp$isocode=="IND",c("year", "gfr", "pab", "gni", "gini", "regSA", "regSSA")], by=c("year"))

## Append the two datasets
dp <- rbind(dp, di[,names(dp)])
dim(dp)



#### check density of variables in estimation and prediction models ####

load(file="../data/data_missing_2020_05_27.RData") ## Load model data

studies <- mutate(studies, sba=if_else(sba>1, sba/100, sba))

## NGA
dp[which(is.na(dp$sba)),]
dp <- mutate(dp, sba=if_else(sba>1, sba/100, sba))
dp[which(is.na(dp$sba)),]


#q()

fcd <- function(i) ggplot() + geom_density(studies, mapping=aes(x=get(vpr[i]))) + geom_density(dp, mapping=aes(x=get(vpr[i])), color="red") + labs(title=paste("Density of",vpr[i],"in both data sets"), x=paste(vpr[i]))
#fcd(13)


#save(dp, file="../results/data_prediction_test.RData")

summary(dp$sba)
dp[which(is.na(dp$sba)),]


#################################################X
####
####      PREDICTIONS with VA models   ##########X


# Load VA model
load("../results/VA_225_0.14.RData") # models
va_225_0.14$output$BUGSoutput$sims.list <- NULL
va_225_0.14$output$BUGSoutput$sims.matrix <- NULL

# matrices of coeficients
mcva <- f.par(va_225_0.14)

# define relevant VA random terms for each study
RVA=sapply(dp$isocode, function(x) c(filter(studies, isocode==x, natrep==T)$reterm))

pva.e <- f.pr2(mcva, PD=mutate(dp, per.early=1, per.late=0, id=paste(isocode,year,sep=".")), PE="early", RT=RVA) # early deaths
pva.l <- f.pr2(mcva, PD=mutate(dp, per.early=0, per.late=1, id=paste(isocode,year,sep=".")), PE="late", RT=RVA) # late deaths

summary(pva.l$PR["NGA.2010","preterm",])
summary(dp %>% filter(isocode=="NGA" & year==2010))

################################################X
####
####      PREDICTIONS with VR models   ##########X

# define relevant VA random terms for each study
RVR=sapply(dp$isocode, function(x) NULL)

# Distribution of early deaths
load("../results/VR1.RData")
vr1_325$output$BUGSoutput$sims.list <- NULL
vr1_325$output$BUGSoutput$sims.matrix <- NULL

mcvr1 <- f.par(vr1_325)
pvr.e <- f.pr2(mcvr1, PD=mutate(dp, id=paste(isocode,year,sep=".")), PE="early", RT=RVR)

#### Predictions of late deaths
load("../results/VR2.RData")
vr2_400$output$BUGSoutput$sims.list <- NULL
vr2_400$output$BUGSoutput$sims.matrix <- NULL

mcvr2 <- f.par(vr2_400)
pvr.l <- f.pr2(mcvr2, PD=mutate(dp, id=paste(isocode,year,sep=".")), PE="late", RT=RVR)


## Save predictions
#save(dp, RVA, pva.e, pva.l, pvr.e, pvr.l, file="../results/prediction_distributions2_test.Rdata")

length(pva.l)

###  Calculate CI for predictions
###
#load("../results/prediction_distributions2_test.Rdata")

## combine predictions early :and late in certain proportions
pe <- 0.74
# For VA model
pva.a <- list(Model="VA", Period=paste0("All (",pe,"e)"), pe=pe, PF=(pe*pva.e$PF+(1-pe)*pva.e$PF), PR=(pe*pva.e$PR+(1-pe)*pva.e$PR))
# For VR model
pvr.a <- list(Model="VR", Period=paste0("All (",pe,"e)"), pe=pe, PF=(pe*pvr.e$PF+(1-pe)*pvr.e$PF), PR=(pe*pvr.e$PR+(1-pe)*pvr.e$PR))


###  Combine VA and VR predictions
###
str(pva.a$PR)
str(pvr.a$PR)
## For each country/year crate a weight of VA dep
dp$isoyear <- paste(dp$isocode, dp$year, sep=".")
dp$wei.vavr <- sapply(dp$nmr, function(x) max(0, min(1, (x - 10)/10 )))
#plot(dp$nmr, dp$wei.vavr)

# causes in each model
ca <- dimnames(pva.a$PF)[[2]]
cr <- dimnames(pvr.a$PF)[[2]]
cc <- unique(c(cr,ca)) # All causes
cc




## Complete matrices with non-existing causes
# in VA matrices
PFa <- abind(pva.a$PF, array(0, dim=c(dim(pva.a$PF)[1], length(cc[!(cc %in% ca)]), dim(pva.a$PF)[3])), along=2)
dimnames(PFa)[[2]] <- c(ca, cc[!(cc %in% ca)])
PFa <- PFa[dp$isoyear,cc,]
PFa[1:3,,1:2]
PRa <- abind(pva.a$PR, array(0, dim=c(dim(pva.a$PR)[1], length(cc[!(cc %in% ca)]), dim(pva.a$PR)[3])), along=2)
dimnames(PRa)[[2]] <- c(ca, cc[!(cc %in% ca)])
PRa <- PRa[dp$isoyear,cc,]
PRa[1:3,,1:2]
# in VR matrices
PFr <- abind(pvr.a$PF, array(0, dim=c(dim(pvr.a$PF)[1], length(cc[!(cc %in% cr)]), dim(pvr.a$PF)[3])), along=2)
dimnames(PFr)[[2]] <- c(cr, cc[!(cc %in% cr)])
PFr <- PFr[dp$isoyear,cc,]
PFr[1:3,,1:2]
PRr <- abind(pvr.a$PR, array(0, dim=c(dim(pvr.a$PR)[1], length(cc[!(cc %in% cr)]), dim(pvr.a$PR)[3])), along=2)
dimnames(PRr)[[2]] <- c(cr, cc[!(cc %in% cr)])
PRr <- PRr[dp$isoyear,cc,]
PRr[1:3,,1:2]

## array of weights
WM <- array(dp$wei.vavr, dim=dim(PFa))
WM[123:125,,1:2]

## Combination
pco.a <- list(Model="VAVR", Period=pva.a$Period, pe=pva.a$pe, PF=(PFa*WM + PFr*(1-WM)), PR=(PRa*WM + PRr*(1-WM)))


summary(pco.a$PR)
#dim(pco.a$PR)



# Calculate CI for predictions
pci.a <- rbind(f.pci2(pva.a), f.pci2(pvr.a), f.pci2(pco.a))

table(pci.a$cause, pci.a$model)

table(pci.a$model, pci.a$type, useNA="ifany")

dim(pci.a)

names(pci.a)
head(pci.a %>% filter(iso=="NGA" & cause=="preterm" & type=="random"), n=30)

dim(pva.a$PR)
dim(pva.a$PR)
head(pva.a$PR["NGA.2011",,1])
summary(pva.a$PR["NGA.2010","preterm",])




write.dta(pci.a, file=paste0("../results/PCIarc2.dta"))

