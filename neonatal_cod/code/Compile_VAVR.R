
rm(list=ls())
library(foreign)
library(readstata13)
library(reshape)

getwd()
label<- "neonatal"

#After model averaging
#causes are stored as fractions
estL<- read.dta13("../results/PCIarc2.dta")
estL<- estL[which(estL$iso !="CHN"),]
estL$iso3<- estL$iso
table(estL$model)
summary(estL$me[which(estL$iso=="IBih")])
estL<- estL[which(estL$year>=2000 & estL$type=="random" & estL$model=="VAVR"),
            c("iso3","year","cause","me")]

#remove indian states
table(nchar(estL$iso3))
estL<- estL[which(nchar(estL$iso3) ==3),]
length(unique(estL$iso3))
estL$iso.year<- paste0(estL$iso3, ".", estL$year)
estL$iso3<-estL$year<-NULL
est<-reshape(estL, direction="wide", timevar="cause", idvar="iso.year")
est<- as.data.frame(est)
est$iso3<- substr(est$iso.year, 1,3)
est$year<- substr(est$iso.year, 5,8)
names(est)<- gsub("me.","",names(est))



#neonatal VA causes
#intrapartum congenital preterm sepsis penumonia diarrhoea tetanus other

#neonatal VR causes
#intrapartum congenital preterm sepsis pneumonia injuries other


#summary(est[which(est$other==0),c("othergroup1","otherncd")])
#est$other<- ifelse(est$other==0, est$othergroup1+est$otherncd, est$other)
#est$othergroup1<- est$otherncd<-NULL



names(est)
vdt<- names(est)[c(2:10)]
vdt

table(rowSums(est[,vdt]))
summary(rowSums(est[,vdt]))

head(est)
is.character(est$iso3)



for(cause in vdt) {
  est[,paste0("f_",cause)]<- est[,cause]
}
f_vdt<- paste0("f_",vdt)
summary(est)
table(rowSums(est[, f_vdt]))

#Envelopes
env<- read.dta13("../data/envelopes_IGME_20200820.dta")
dim(env)
head(env)
summary(env)
names(env)




temp2<- merge(est, env, by=c("iso3","year"))
summary(temp2[,f_vdt])
unique(temp2$iso3)

#remove PSE
est<- temp2[which(temp2$iso3 !="PSE"),]

##########################################################################

summary(est$nmr)
names(est)[order(names(est))]
est$model<- #ifelse(est$iso3 %in% gvr$iso3, "Good VR",
            ifelse(est$nmr>=20, "VA", ifelse(est$nmr<=10, "VR","VA/VR mix"))
table(est$model)

write.dta(est, file=paste0("../results/Estimates_",label,".dta"))
summary(est[,f_vdt])
q("no")
