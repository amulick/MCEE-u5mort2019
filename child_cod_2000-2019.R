#This "master file" produces WHO-MCEE country-level child cause of death estimates for 2000-2019.

#*See "Progress towards the Sustainable Development Goals:  Global, regional,
#      and national causes of under-five mortality during 2000-2019" for a description of methods used to produce these estimates.

#########################################################################################################
#  Cause cateogries output by this program
# 	CH2		HIV/AIDS
# 	CH3		Diarrhoeal diseases
# 	CH5		Tetanus
# 	CH6		Measles
# 	CH7		Meningitis/encephalitis
# 	CH8		Malaria
# 	CH9		Acute respiratory infections
# 	CH10	Prematurity
# 	CH11	Birth asphyxia and birth trauma (also referred to as intrapartum-related complications)
# 	CH12	Sepsis and other infectious conditions of the newborn
# 	CH13	Other Group 1 (1-59 months only)
# 	CH15	Congenital anomalies
# 	CH16	Other noncommunicable diseases (1-59 months only)
# 	CH17	Injuries
#   CH18 - TB
#   CH19 - Other Group 1 and Other noncommunicable (neonatal and under five only)
#########################################################################################################

rm(list=ls())
library(foreign)
library(readstata13)
library(dbplyr)
library(tidyverse)
library(dplyr)
library(magrittr)
options(width=60)


##########################################################################################################
##########################################################################################################
##  A. PREPARE ESTIMATES FOR CAUSES OF NEONATAL DEATHS
##########################################################################################################
##########################################################################################################

#Combine neonatal inputs
#Cause fractions from countries with high quality registration
#n.goodvr<-read.dta13("../data/neonatal/good_vr_neonates_20201012.dta")
n.goodvr<- read.dta13("neonatal_cod/data/good_vr_neonates_20210106.dta")
n.goodvr$iso3<- n.goodvr$isocode
n.goodvr$model<- "Good VR"


#n.vavr<-read.dta13( "../data/neonatal/vavr_neonatal_20200813.dta")
n.vavr<- read.dta13("neonatal_cod/results/Estimates_neonatal.dta")
n.vavr<- n.vavr[which(!n.vavr$iso3 %in% n.goodvr$iso3),]


n.vavr<- n.vavr[,c("iso3","year","model","f_preterm",     "f_intrapartum", "f_injuries",
        "f_congenital",  "f_sepsis", "f_pneumonia",   "f_diarrhoea",  "f_tetanus",     "f_other")]
n.dat<- merge(n.goodvr, n.vavr, by=c("iso3","year","model"), all=TRUE)
n.dat[which(n.dat$iso3=="MCO"),c("iso3","year","f_preterm","f_sepsis","f_other","model")]

#Merge envelopes and program estimates
other_inputs<- read.dta("other_inputs/other_inputs_20210203.dta")
n.dat<- merge(n.dat, other_inputs, by=c("iso3" ,"year"))
n.dat<-n.dat %>% filter(iso3 !="CHN")


#*Map to numbered variables for convenience
n.dat <- n.dat %>% mutate(frac3  = f_diarrhoea,
frac5  = f_tetanus,
frac9  = f_pneumonia,
frac10 = f_preterm,
frac11 = f_intrapartum,
frac12 = f_sepsis,
frac15 = f_congenital,
frac17 = f_injuries,
frac_oth=f_other,
frac3 = ifelse(model=="VR",0, frac3))
n.dat<- as.data.frame(n.dat)

#Check
tot<-rowSums(n.dat[,c("frac3","frac5","frac9","frac10","frac11","frac12","frac15","frac17","frac_oth")])
temp<- rowSums(n.dat[,c("frac3","frac5","frac9","frac10","frac11","frac12","frac15","frac17","frac_oth")])
tapply(temp, n.dat$model, FUN=summary)




#Crisis neonatal deaths attributed to all causes in Venezuela health crisis
n.dat$nnd<- ifelse(n.dat$iso3=="VEN",
                  n.dat$igme_nnd_pcrisis - n.dat$hivn,  n.dat$nnd)
#tapply(n.dat$nnt/n.dat$nnd, n.dat$model, FUN=summary)
#tapply(n.dat$frac5, n.dat$model, FUN=summary)
#tapply(n.dat$frac5*n.dat$nnd, n.dat$model, FUN=summary)



#Map causes for modeled countries
n.dat <- n.dat %>% mutate(neo2  = ifelse(model=="Good VR", hiv_filled*(nnd+hivn), hivn),
                          neo3  = ifelse(model=="Good VR", 0, frac3*(nnd+hivn)),
                          neo5  = case_when(model=="Good VR"~ tetanus_filled*(nnd+hivn),
                                            model=="VA" ~ frac5*(nnd+hivn),
                                            TRUE ~ frac5*(nnd+hivn)),
                          neo6  = 0,
                          neo7  = ifelse(model=="Good VR", meningitis_filled*(nnd+hivn),
                                         men_sepsis_frac*frac12*(nnd+hivn)),
                          neo8  = ifelse(model=="Good VR", malaria_filled*(nnd+hivn),0),
                          neo9  = ifelse(model=="Good VR", (pneumonia_filled)*(nnd+hivn), frac9*(nnd+hivn)),
                          neo10  = ifelse(model=="Good VR", preterm_filled*(nnd+hivn), frac10*(nnd+hivn)),
                          neo11  = ifelse(model=="Good VR", intrapartum_filled*(nnd+hivn), frac11*(nnd+hivn)),
                          neo12  = case_when(model=="Good VR"~ sepsis_filled*(nnd+hivn),
                                          model=="VA"~ (1-men_sepsis_frac)*frac12*(nnd+hivn) -hivn/2,
                                          TRUE ~ (1-men_sepsis_frac)*frac12*(nnd+hivn) -hivn/2),
                          neo15  = ifelse(model=="Good VR", congenital_filled*(nnd+hivn), frac15*(nnd+hivn)))


n.dat<- n.dat %>% mutate(neo19  = case_when(model=="Good VR"~ (other1_filled+other2_filled)*(nnd+hivn) ,
                                            model =="VA" ~ frac_oth*(1-vamcm_otherneo17)*(nnd+hivn) -hivn/2,
                                            TRUE~ frac_oth*(nnd+hivn) -hivn/2),
                         neo17  = case_when(model=="Good VR"~ injuries_filled*(nnd+hivn),
                                           model=="VA" ~ frac_oth*vamcm_otherneo17*(nnd+hivn) + frac17*(nnd+hivn) ,
                                           TRUE~ frac17*(nnd+hivn)),
                         neo18 = 0) #TB




#Use UNAIDS estimates for HIV in Jamaica
n.dat$neo2<- ifelse(n.dat$iso3=="JAM",n.dat$hivn, n.dat$neo2)
neotot <- rowSums(n.dat[, paste0("neo",c(2:3,5:12,19,17:18))])


#Check sums against envelope
summary(n.dat[, paste0("neo", c(2:3,5:12,15,19,17:18)) ])
neotot <- rowSums(n.dat[, paste0("neo",c(2:3,5:12,15,19,17:18))])



tapply(neotot-(n.dat$nnd+n.dat$hivn), n.dat$model, FUN=summary)
n.dat$diff<-neotot-(n.dat$nnd+n.dat$hivn)
tapply((n.dat$nnd+n.dat$hivn), n.dat$model, FUN=summary)
table(n.dat$model)
summary(n.dat$diff/n.dat$nnd)
n.dat[which(abs(n.dat$diff) > 100),c("iso3","year","vamcm_otherneo13",
                                     "vamcm_otherneo16","vamcm_otherneo17", "frac17","nnd")]
n.dat$neometh<- n.dat$model


for (X in c(2:3, 5:12,15,19,17:18)) {
    print(X)
  n.dat[, paste0("neo",X)]<-        ifelse(neotot==0 | (n.dat$nnd+n.dat$hivn==0) , n.dat[, paste0("neo",X)],
                                           n.dat[, paste0("neo",X)] /neotot*(n.dat$nnd+n.dat$hivn) )
}

#Check envelopes
neotot <- rowSums(n.dat[, paste0("neo",c(2:3,5:12,15, 19, 17:18))])
tapply(neotot/(n.dat$igme_nnd_pcrisis), n.dat$model, FUN=summary)
tapply(neotot - (n.dat$igme_nnd), n.dat$model, FUN=summary)




####################################################################################################
#For countries with modeled estimates, using tetanus program numbers
for(X in paste0("neo",c(3,7,9:12,15,19,17:18))) {
  #replace postX = postX * (pnd-mal-post6)/(pnd-post8-post6) if malwho==1 & pnmeth=="vamcm"
  n.dat[,X]<- ifelse(n.dat$model=="VR" , n.dat[,X]*(n.dat$nnd - n.dat$nnt) /
                       (n.dat$nnd - n.dat$neo5),
                     n.dat[,X])
}
n.dat$neo5<- ifelse(n.dat$model=="VR",n.dat$nnt, n.dat$neo5)

neotot <- rowSums(n.dat[, paste0("neo",c(2:3,5:12,15,19,17:18))])
tapply(neotot/(n.dat$igme_nnd_pcrisis), n.dat$model, FUN=summary)
tapply(neotot - (n.dat$igme_nnd), n.dat$model, FUN=summary)
n.dat$diff<- neotot - (n.dat$igme_nnd)
n.dat[which(abs(n.dat$diff)> .1),c("iso3","year","neo5","nnt","frac5","nnd","model","diff")]

#Include crisis deaths for Zika as congenital
n.dat$neo15<- ifelse(n.dat$iso3=="BRA" & n.dat$year==2016, n.dat$neo15+ n.dat$igme_nncrisis, n.dat$neo15)


#Yemen
#Include crisis deaths for Ebola as other infections/malaria
#n.dat[which(n.dat$iso3=="LBR"),c("year","igme_nnd","igme_nncrisis")]
#n.dat[which(n.dat$iso3=="GIN"),c("year","igme_nnd","igme_nncrisis")]
#n.dat[which(n.dat$iso3=="SLE"),c("year","igme_nnd","igme_nncrisis")]
#n.dat[which(n.dat$iso3=="YEM"),c("year","igme_nnd","igme_nncrisis")]

#Ebola
n.dat$neo19 <- ifelse( n.dat$iso3 %in% c("SLE","LBR","GIN"),
                       n.dat$neo19+ n.dat$igme_nncrisis, n.dat$neo19)

#Include other crisis deaths as injuries
#For Yemen include as fraction of other group 1
n.dat$neo17<- ifelse( !(n.dat$iso3=="BRA" & n.dat$year==2016) & !n.dat$iso3 %in% c("VEN","SLE","LBR","GIN"),
                      n.dat$neo17+ n.dat$igme_nncrisis*(1-n.dat$frac_pncrisis_maln), n.dat$neo17)
n.dat$neo19<- ifelse( !(n.dat$iso3=="BRA" & n.dat$year==2016) & !n.dat$iso3 %in% c("VEN","SLE","LBR","GIN"),
                      n.dat$neo19+ n.dat$igme_nncrisis*(n.dat$frac_pncrisis_maln), n.dat$neo19)


neotot <- rowSums(n.dat[, paste0("neo",c(2:3,5:12,15,19,17:18))])
tapply(neotot/(n.dat$igme_nnd_pcrisis), n.dat$model, FUN=summary)
tapply(neotot - (n.dat$igme_nnd_pcrisis), n.dat$model, FUN=summary)

#temp<-neotot - (n.dat$igme_nnd_pcrisis)
#n.dat[which(temp< -.1),]


#include China
chn<-read.dta13("neonatal_cod/data/china_neonatal_extrapolated_20200824.dta")
chn$neometh<- "Other"
chn<- merge(chn, other_inputs[,c("iso3","year","hivn")], by=c("iso3","year"))
dim(chn)
names(chn)
chn[,c("year","nnd","hivn")]
chn$neo19 <- chn$neo13+chn$neo16


#Update HIV
for(X in paste0("neo",c(3,5,7,9:11,12,15,17,19))) {
  #replace postX = postX * (pnd-mal-post6)/(pnd-post8-post6) if malwho==1 & pnmeth=="vamcm"
  chn[,X]<-  chn[,X]*(chn$nnd - chn$hivn) / (chn$nnd - chn$neo2)

}
chn$neo2<- chn$hivn


chn$neo18<-0
chn$igme_nnd_pcrisis<- chn$nnd
n.dat$whoname<- n.dat$country
var<-c("iso3","whoname","year", paste0("neo",c(2:3,5:12,15,19,17:18)),"nnd","neometh","igme_nnd_pcrisis")
n.dat<- rbind(n.dat[,var],chn[,var])

summary(n.dat[,paste0("neo",c(2:3,5:12,15,19,17:18))])
summary(rowSums(n.dat[,paste0("neo",c(2:3,5:12,15,19,17:18))]))
summary(rowSums(n.dat[,paste0("neo",c(2:3,5:12,15,19,17:18))]) - n.dat$nnd)
n.dat$nnd<-rowSums(n.dat[,paste0("neo",c(2:3,5:12,15,19,17:18))])

#Check against IGME envelopes
tapply(n.dat$nnd, list(n.dat$year), FUN=sum)
tapply(n.dat$igme_nnd_pcrisis, list(n.dat$year), FUN=sum)

#n.dat[which(n.dat$igme_nnd_pcrisis > n.dat$nnd), c("iso3","year","nnd","igme_nnd_pcrisis")]
summary(n.dat$igme_nnd_pcrisis - n.dat$nnd)
summary(n.dat$nnd - rowSums(n.dat[,paste0("neo",c(2:3,5:12,15,19,17:18))]) )
summary(n.dat[,paste0("neo",c(2:3,5:12,15,19,17:18))])






####################################################################################################
####################################################################################################
##  B. PREPARE ESTIMATES FOR CAUSES OF POSTNEONATAL DEATHS
####################################################################################################
####################################################################################################

#pn.goodvr<- read.dta13("../data/postneonatal/postneonatal_goodvr_1990-2019_20200515.dta")
pn.goodvr<- read.dta13("postneonatal_cod/data/postneonatal_goodvr_1990_20200915.dta")
pn.goodvr<- pn.goodvr %>% filter(year>=2000)
pn.goodvr$pnmeth<-"Good VR"

#results VA and VR models, after model averaging, saved in "est"
load("postneonatal_cod/results/Model_averaging.RData")
#est<- read.dta("postneonatal_cod/results/Prediction2_model_averaging_JP.dta")

pn.vavr<- as.data.frame(est)
pn.vavr$iso3<- substr(rownames(pn.vavr),1,3)
pn.vavr$year<- substr(rownames(pn.vavr),5,8)
summary(rowSums(pn.vavr[,c("diarrhea","pneumonia","other","congenital","neonatal",
                           "meningitis","othergroup1","otherncd","injuries","malaria")]))

#All countries except China
p.dat<- merge(pn.goodvr, pn.vavr, by=c("iso3","year"), all=TRUE)
#n.dat$iso3[ !n.dat$iso3 %in% p.dat$iso3]

p.dat<- merge(p.dat, other_inputs, by=c("iso3" ,"year"))

p.dat$pnmeth<- ifelse(is.na(p.dat$pnmeth),
                      ifelse(p.dat$u5mr<=25, "VR",
                             ifelse(p.dat$u5mr>=35, "VA","VA/VR"
                             )), p.dat$pnmeth)


p.dat$malaria<- ifelse(is.na(p.dat$malaria),0, p.dat$malaria)
p.dat$malwho<- ifelse(is.na(p.dat$malwho),0, p.dat$malwho)


#Crisis deaths attributed to all causes in Venezuela health crisis
p.dat$pnd<- ifelse(p.dat$iso3=="VEN",
                   p.dat$igme_pnd_pcrisis - p.dat$hivp,  p.dat$pnd)


#check envelopes
tapply(rowSums(other_inputs[,c("pnd","hivp")]), other_inputs$year, FUN=sum)
tapply(other_inputs[,c("igme_pnd")], other_inputs$year, FUN=sum)



p.dat <- p.dat %>% mutate(post2  = ifelse(pnmeth=="Good VR", post2*(pnd+hivp), hivp),
                          post3  = case_when(pnmeth=="Good VR"~ post3*(pnd+hivp),
                                             pnmeth=="VR"~ diarrhea*(pnd-measin),
                                             pnmeth %in% c("VA", "VA/VR")~ diarrhea*(pnd-measin)),
                          post5  = ifelse(pnmeth=="Good VR", post5*(pnd+hivp), pnt),
                          post6  = ifelse(pnmeth=="Good VR", post6*(pnd+hivp), measin),
                          post7  = case_when(pnmeth=="Good VR"~ post7*(pnd+hivp),
                                          pnmeth=="VR"~ meningitis*(pnd-measin),
                                          TRUE ~ meningitis*(pnd-measin)),
                          post8  = ifelse(pnmeth=="Good VR", post8*(pnd+hivp),
                                          malaria*(pnd-measin) ),
                          post9  = case_when(pnmeth=="Good VR"~ (post9+post4)*(pnd+hivp),
                                             pnmeth=="VR"~pneumonia*(pnd-measin),
                                             TRUE~pneumonia*(pnd-measin)))



p.dat<- p.dat %>% mutate( post10  = case_when(pnmeth=="Good VR"~ post10*(pnd+hivp),
                                              pnmeth=="VR"~ post_peri_prem_frac*neonatal*(pnd-measin),
                          TRUE~ post_peri_prem_frac*neonatal*(pnd-measin)),
                          post11  = case_when(pnmeth=="Good VR"~ post11*(pnd+hivp),
                                              pnmeth=="VR"~ (1-post_peri_prem_frac)*neonatal*(pnd-measin),
                                                      TRUE~(1-post_peri_prem_frac)*neonatal*(pnd-measin)),
                          post12  = ifelse(pnmeth=="Good VR", post12*(pnd+hivp), 0),
                          post13  = case_when(pnmeth=="Good VR"~ post13*(pnd+hivp) ,
                                              other ==0 ~ (othergroup1 - post8/(pnd-measin))*(pnd-measin)-pnt ,
                                              TRUE~ (post_other_gr1_frac*(other)*(pnd-measin) + othergroup1*(pnd-measin) - pnt )),
                          post15  = case_when(pnmeth=="Good VR"~ post15*(pnd+hivp),
                                           pnmeth=="VR"~congenital*(pnd-measin),
                                           TRUE~congenital*(pnd-measin)),
                          post16  = case_when(pnmeth=="Good VR"~ post16*(pnd+hivp),
                                              other==0 ~ otherncd*(pnd-measin),
                                      TRUE~  (1-post_other_gr1_frac)*(other)* (pnd-measin) + otherncd*(pnd-measin)),
                          post17  = case_when(pnmeth=="Good VR"~ post17*(pnd+hivp),
                                              pnmeth=="VR"~injuries*(pnd-measin),
                                               TRUE~ injuries*(pnd-measin))
                            )


#CHECK ENVELOPES
#doesn't exactly match IGME crisis free because of VEN
posttot <- rowSums(p.dat[, paste0("post",c(2:3,5:13,15:17))])
tapply(posttot, p.dat$year, FUN=sum)
tapply(p.dat$pnd + p.dat$hivp, p.dat$year, FUN=sum)
summary(posttot - p.dat$pnd - p.dat$hivp)

#Cap TB at other group 1
p.dat$tb <- p.dat$tb_resp+p.dat$tb_nonresp



#Residual From countries with good VR
#> summary(est.gvr$post13/(est.gvr$pnd))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.02937 0.04286 0.05275 0.06955 0.39474
residual_post13_f<- 0.05275
p.dat$residual_post13<- p.dat$pnd*residual_post13_f



#> summary(est.gvr$post9/(est.gvr$pnd))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.05386 0.08983 0.11973 0.16598 0.57056
residual_ari_f <- 0.11973
p.dat$residual_ari<- p.dat$pnd*residual_ari_f


#Instances where TB is restricted by LRI or by other group 1 fractions
table((p.dat$tb_resp > p.dat$post9))
table((p.dat$tb_resp > p.dat$post9 - residual_ari_f*p.dat$pnd))
table((p.dat$tb_nonresp > p.dat$post13))
table((p.dat$tb_nonresp > p.dat$post13 - residual_post13_f*p.dat$pnd))




p.dat<- p.dat %>% mutate( tb_nonresp = ifelse(is.na(tb_nonresp),0, tb_nonresp),
                          tb_resp = ifelse(is.na(tb_resp),0, tb_resp))


#squeezing TB and residual ARI
# Non respiratory TB, (out of other CD)
p.dat<- p.dat %>% mutate( tb_nonresp_corr = ifelse(tb_nonresp + residual_post13 > post13,
                                              tb_nonresp/(tb_nonresp + residual_ari)*post13, tb_nonresp),
                          tb_resp_corr = ifelse(tb_resp + residual_ari > post9,
                                           tb_resp/(tb_resp + residual_ari)*post9, tb_resp))



#Add excess from non respiratory TB to respiratory TB
p.dat$tb_resp_plus<- p.dat$tb_resp + (p.dat$tb_nonresp - p.dat$tb_nonresp_corr)
p.dat<- p.dat %>% mutate( tb_resp_plus_corr = ifelse(tb_resp_plus+ residual_ari > post9,
                                           tb_resp_plus/(tb_resp_plus+ residual_ari)*post9, tb_resp_plus))





# Respiratory TB, (out of LRI)
p.dat<- p.dat %>% mutate( residual_ari_corr = ifelse(tb_resp_plus+ residual_ari > post9,
                                              residual_ari/(tb_resp_plus+ residual_ari)*post9, post9 - tb_resp),
                          residual_post13_corr = ifelse(tb_nonresp + residual_post13 > post13,
                                           residual_post13/(tb_nonresp + residual_post13)*post13, post13 - tb_nonresp))


#What fraction is respiratory TB out of all LRI mortality?
#Before addition from non respiratory
summary(p.dat$tb_resp_corr/(p.dat$post9))
tapply(p.dat$tb_resp_corr/(p.dat$post9), p.dat$pnmeth, FUN=summary)

#After addition from non respiratory
summary(p.dat$tb_resp_plus_corr/(p.dat$post9))
tapply(p.dat$tb_resp_plus_corr/(p.dat$post9), p.dat$pnmeth, FUN=summary)


#for a few country years where TB missing (due to boundary changes)
p.dat$post13.org<- p.dat$post13
p.dat$post9.org<- p.dat$post9
p.dat<- p.dat %>% mutate( post13 = post13 - tb_nonresp_corr,
                          post9 = post9 - tb_resp_plus_corr,
                          post18 = tb_resp_plus_corr+tb_nonresp_corr)

#CHECK ENVELOPES
#doesn't exactly match IGME crisis free because of VEN
posttot <- rowSums(p.dat[, paste0("post",c(2:3,5:13,15:18))])
tapply(posttot, p.dat$year, FUN=sum)
tapply(p.dat$pnd + p.dat$hivp, p.dat$year, FUN=sum)
summary(posttot - p.dat$pnd - p.dat$hivp)


tapply(posttot - (p.dat$pnd+p.dat$hivp), p.dat$pnmeth, FUN=summary)
p.dat$posttot<- posttot
p.dat$diff<-p.dat$posttot - (p.dat$pnd+p.dat$hivp)
summary(p.dat$diff)
table(abs(p.dat$diff) >.1, p.dat$malwho, useNA="always")
p.dat[order(-abs(p.dat$diff))[1:20], c("iso3","year","posttot","pnd","igme_pnd",
                                       "hivp","pnmeth","diff","other","othergroup1","post8")]

posttot <- rowSums(p.dat[, paste0("post",c(2:3,5:13,15:18))])
summary(p.dat[, paste0("post",c(2:3,5:13,15:18))])


####################################################################################################
#For countries with modeled estimates, using malaria program numbers
table(p.dat$pnmeth, p.dat$malwho)
for(X in paste0("post",c(3,5,7,9:11,13,15:17))) {
  p.dat[,X]<- ifelse(p.dat$malwho==1 & p.dat$pnmeth !="Good VR", p.dat[,X]*(p.dat$pnd - p.dat$mal - p.dat$post6 - p.dat$post18) /
                       (p.dat$pnd - p.dat$post8 - p.dat$post6 - p.dat$post18),
                        p.dat[,X])
}
tapply(p.dat$post8 - p.dat$mal, p.dat$pnmeth, FUN=summary)
p.dat$post8<- ifelse(p.dat$malwho==1 & p.dat$pnmeth !="Good VR",p.dat$mal, p.dat$post8)
tapply((p.dat$post8 - p.dat$mal)/p.dat$pnd, p.dat$pnmeth, FUN=summary)
tapply(p.dat$post8 - p.dat$mal, p.dat$pnmeth, FUN=function(x) sum(is.na(x)))


#For countries with good VR using malaria program estimates
for(X in paste0("post",c(2,3,5:7,9:11,13,15:18))) {
  p.dat[,X]<- ifelse(p.dat$malwho==1 & p.dat$pnmeth =="Good VR", p.dat[,X]*(p.dat$pnd - p.dat$mal+ p.dat$hivp) /
                       (p.dat$pnd - p.dat$post8 + p.dat$hivp),
                     p.dat[,X])
}
p.dat$post8<- ifelse(p.dat$malwho==1 & p.dat$pnmeth =="Good VR", p.dat$mal, p.dat$post8)

#CHECK ENVELOPES
posttot <- rowSums(p.dat[, paste0("post",c(2:3,5:13,15:18))])
tapply(posttot, p.dat$year, FUN=sum)
tapply(p.dat$pnd+p.dat$hivp, p.dat$year, FUN=sum)
tapply(posttot - (p.dat$pnd+p.dat$hivp), p.dat$pnmeth, FUN=summary)
tapply(posttot - (p.dat$pnd+p.dat$hivp), p.dat$malwho, FUN=summary)
p.dat$posttot<- posttot
p.dat$diff<-p.dat$posttot - (p.dat$pnd+p.dat$hivp)
summary(p.dat$diff)
table(abs(p.dat$diff) >.1, p.dat$malwho, useNA="always")
p.dat[order(-abs(p.dat$diff))[1:20], c("iso3","year","posttot","pnd","igme_pnd",
                                       "hivp","pnmeth","diff","post8","mal","post2")]

posttot <- rowSums(p.dat[, paste0("post",c(2:3,5:13,15:18))])
summary(p.dat[, paste0("post",c(2:3,5:13,15:18))])



#summary(p.dat$post8[which(p.dat$iso3=="TUV")])
summary(p.dat[, paste0("post",c(2:3,5:13,15:18))])
p.dat$post4<- NULL


#Use UNAIDS estimates for HIV in Jamaica
for(X in paste0("post",c(3,7,9:13,15:18))) {
  p.dat[,X]<- ifelse(p.dat$iso3=="JAM" , p.dat[,X]*(p.dat$igme_pnd - p.dat$hivp) /
                       (p.dat$igme_pnd - p.dat$post2),
                     p.dat[,X])
}
p.dat$post2<- ifelse(p.dat$iso3=="JAM",p.dat$hivp, p.dat$post2)
p.dat[which(p.dat$iso3=="JAM"),c("post2","hivp","pnd","igme_pnd")]


#CHECK ENVELOPES
posttot <- rowSums(p.dat[, paste0("post",c(2:3,5:13,15:18))])
tapply(posttot, p.dat$year, FUN=sum)
tapply(p.dat$pnd+p.dat$hivp, p.dat$year, FUN=sum)
tapply(p.dat$igme_pnd, p.dat$year, FUN=sum)

tapply(posttot - (p.dat$pnd+p.dat$hivp), p.dat$pnmeth, FUN=summary)
p.dat$posttot<- posttot
p.dat$diff<-p.dat$posttot - (p.dat$pnd+p.dat$hivp)
summary(p.dat$diff)
table(abs(p.dat$diff) >.1, p.dat$malwho, useNA="always")
p.dat[order(-abs(p.dat$diff))[1:20], c("iso3","year","posttot","pnd","igme_pnd","mal","othergroup1","other",
                                       "post8","malwho","pnmeth","diff","measin","post6","u5mr","igme_pncrisis","post8")]




#China
ch<- read.dta("postneonatal_cod/data/china_postneonatal_extrapolated_20200824.dta")
for(X in c(3, 5:13, 15:17)) {
  ch[, paste0("post",X)]<- ch[, paste0("post",X)]/ch$pnd
}
ch$pnd<- NULL
ch$pnmeth="Empirical"

chn<- merge(ch, other_inputs, by=c("iso3","year"))
for(X in c(3, 5:13, 15:17)) {
  chn[, paste0("post",X)]<- chn[, paste0("post",X)]*(chn$pnd+chn$hivp)
}

for(X in paste0("post",c(3,5,7,9:11,13,15:17))) {
  chn[,X]<- ifelse(chn$malwho==1, chn[,X]*(chn$pnd +chn$hivp - chn$mal - chn$measin - chn$hivp) /
                       (chn$pnd +chn$hivp- chn$post8 - chn$post6 - chn$post2),
                        chn[,X])
}
chn$post8<- ifelse(chn$malwho==1, chn$mal, chn$post8)
chn$post2<- chn$hivp
chn$post6<- chn$measin



#Include TB
chn$residual_post13<- chn$pnd*residual_post13_f
chn$residual_ari<- chn$pnd*residual_ari_f


#squeezing TB and residual ARI
chn<- chn %>% mutate( tb_nonresp_corr = ifelse(tb_nonresp + residual_post13 > post13,
                                              tb_nonresp/(tb_nonresp + residual_ari)*post13, tb_nonresp),
                          tb_resp_corr = ifelse(tb_resp + residual_ari > post9,
                                           tb_resp/(tb_resp + residual_ari)*post9, tb_resp))

chn<- chn %>% mutate( residual_ari_corr = ifelse(tb_resp + residual_ari > post9,
                                              residual_ari/(tb_resp + residual_ari)*post9, post9 - tb_resp),
                          residual_post13_corr = ifelse(tb_nonresp + residual_post13 > post13,
                                           residual_post13/(tb_nonresp + residual_post13)*post13, post13 - tb_nonresp))

#Add excess from non respiratory TB to respiratory TB
chn$tb_resp_plus<- chn$tb_resp + (chn$tb_nonresp - chn$tb_nonresp_corr)
chn<- chn %>% mutate( tb_resp_plus_corr = ifelse(tb_resp_plus+ residual_ari > post9,
                                           tb_resp_plus/(tb_resp_plus+ residual_ari)*post9, tb_resp_plus))
chn<- chn %>% mutate( post13 = post13 - tb_nonresp_corr,
                          post9 = post9 - tb_resp_plus_corr,
                          post18 = tb_resp_plus_corr+tb_nonresp_corr)


var<- c("iso3","year","whoreg6",paste0("post",c(2:3, 5:13, 15:18)), "crisis",
        "measin","measout","pnmeth","lb","u5mr","igme_pnd_pcrisis",
        "frac_pncrisis_maln","igme_pncrisis", "igme_pnd","hivp","pnd","malwho","mal",
        "tb_resp","tb_nonresp","tb_resp_corr","tb_nonresp_corr","pfpr","tb_resp_plus","tb_resp_plus_corr",
        "igme_u5mr_pcrisis")
p.dat$year<- as.numeric(as.character(p.dat$year))
temp<- left_join(p.dat, other_inputs[,c("iso3","year")], by=c("iso3","year"))
p.dat<- temp
for(v in var){
    print(v)
    summary(p.dat[,v])
    cat("ch  ",v,"\n")
        summary(chn[,v])
}

p.dat<- rbind(p.dat[,var], chn[,var])
names(p.dat)
summary(p.dat$post8)


#Check against IGME envelopes
summary(p.dat[,paste0("post",c(2:3,5:13,15:18))])
summary(rowSums(p.dat[,paste0("post",c(2:3,5:13,15:18))]))
summary(rowSums(p.dat[,paste0("post",c(2:3,5:13,15:18))]) - (p.dat$pnd+p.dat$hivp))

p.dat$posttot<-rowSums(p.dat[,paste0("post",c(2:3,5:13,15:18))])
p.dat$diff<-p.dat$posttot - (p.dat$pnd+p.dat$hivp)
p.dat$diff<-p.dat$posttot - (p.dat$igme_pnd)
summary(p.dat$diff)
table(abs(p.dat$diff) >.1, p.dat$malwho, useNA="always")
names(p.dat)
p.dat[order(-abs(p.dat$diff))[1:20], c("iso3","year","posttot","pnd","igme_pnd","mal","diff")] #,#"othergroup1","other",
                                       #"post8","malwho","pnmeth","diff","measin","post6","u5mr","igme_pncrisis","post8")]

summary((p.dat$pnd+p.dat$hivp) - (p.dat$igme_pnd))
#p.dat$pnd<-rowSums(p.dat[,paste0("post",c(2:3,5:13,15:18))])


#CHECK ENVELOPES
tapply(p.dat$measout, p.dat$year, FUN=sum)
summary(p.dat[,paste0("post",c(2:3,5:13,15:18))])
posttot <- rowSums(p.dat[, paste0("post",c(2:3,5:13,15:18))])
tapply(posttot, p.dat$year, FUN=sum)
tapply(p.dat$igme_pnd, p.dat$year, FUN=sum)
tapply(p.dat$pnd+p.dat$hivp, p.dat$year, FUN=sum)

tapply(posttot - (p.dat$pnd+p.dat$hivp), p.dat$pnmeth, FUN=summary)
p.dat$posttot<- posttot
p.dat$diff<-p.dat$posttot - (p.dat$pnd+p.dat$hivp)
summary(p.dat$diff)

#table(abs(p.dat$diff) >.1, p.dat$malwho, useNA="always")
p.dat[order(-abs(p.dat$diff))[1:20], c("iso3","year","posttot","pnd","igme_pnd",
                                       "hivp","pnmeth","diff","post8")]


#include measles outside envelopes
p.dat$post6<- ifelse(p.dat$pnmeth !="Good VR", p.dat$post6+p.dat$measout, p.dat$post6)

#include crisis deaths
summary(p.dat[,paste0("post",c(2:3, 5:13, 15:18))])
names(p.dat)

#CHECK ENVELOPES
tapply(p.dat$measout, p.dat$year, FUN=sum)
summary(p.dat[,paste0("post",c(2:3,5:13,15:18))])
posttot <- rowSums(p.dat[, paste0("post",c(2:3,5:13,15:18))])
tapply(posttot, p.dat$year, FUN=sum)
tapply(p.dat$igme_pnd, p.dat$year, FUN=sum)
tapply(p.dat$pnd+p.dat$hivp + p.dat$measout, p.dat$year, FUN=sum)
tapply(posttot - (p.dat$pnd+p.dat$hivp + p.dat$measout), p.dat$year, FUN=sum)
p.dat$posttot<- posttot
p.dat$diff<-p.dat$posttot - (p.dat$pnd+p.dat$hivp+p.dat$measout)
summary(p.dat$diff)
p.dat[order(-abs(p.dat$diff))[1:20], c("iso3","year","posttot","pnd","igme_pnd",
                                       "measout","diff")]


p.dat$posttot<- posttot
p.dat$diff<-p.dat$posttot - (p.dat$pnd+p.dat$hivp +ifelse(p.dat$pnmeth=="Good VR", 0,p.dat$measout))
summary(p.dat$diff)
p.dat[order(-abs(p.dat$diff))[1:20], c("iso3","year","posttot","pnd","igme_pnd",
                                       "hivp","pnmeth","diff","post8")]

#Include crisis deaths for Zika as congenital
p.dat$post15<- ifelse(p.dat$iso3=="BRA" & p.dat$year==2016, p.dat$post15+ p.dat$igme_pncrisis, p.dat$post15)


#Include crisis deaths for Ebola as other infections/malaria
names(p.dat)
p.dat[which(p.dat$iso=="LBR"),c("year","igme_pnd","igme_pncrisis","measout","measin")]
p.dat[which(p.dat$iso=="GIN"),c("year","igme_pnd","igme_pncrisis","measout","measin")]
p.dat[which(p.dat$iso=="SLE"),c("year","igme_pnd","igme_pncrisis","measout","measin")]

p.dat$post13<- ifelse( p.dat$iso3 %in% c("SLE","LBR","GIN"),
                       p.dat$post13+ p.dat$igme_pncrisis*(0.78), p.dat$post13)
p.dat$post8<- ifelse( p.dat$iso3 %in% c("SLE","LBR","GIN"),
                       p.dat$post8+ p.dat$igme_pncrisis*(0.22), p.dat$post8)


#Include other crisis deaths as injuries
#and partially due to malnutrition in Yemen
p.dat$post17<- ifelse( !(p.dat$iso3=="BRA" & p.dat$year==2016) &
                         !p.dat$iso3 %in% c("VEN","SLE","LBR","GIN"),
                       p.dat$post17+ p.dat$igme_pncrisis*(1-p.dat$frac_pncrisis_maln), p.dat$post17)
p.dat$post13<- ifelse( !(p.dat$iso3=="BRA" & p.dat$year==2016) & !p.dat$iso3 %in% c("VEN","SLE","LBR","GIN"),
                       p.dat$post13+ p.dat$igme_pncrisis*(p.dat$frac_pncrisis_maln), p.dat$post13)

######################################################################
#Check against IGME envelopes
summary(p.dat[,paste0("post",c(2:3,5:13,15:18))])
tapply(rowSums(p.dat[,paste0("post",c(2:3,5:13,15:18))]), p.dat$year, FUN=sum)
#this should match
tapply(p.dat$igme_pnd_pcrisis + ifelse(p.dat$pnmeth=="Good VR", 0,p.dat$measout), p.dat$year, FUN=sum)
#and match
tapply(p.dat$pnd+p.dat$hivp +ifelse(p.dat$pnmeth=="Good VR", 0,p.dat$measout) +
         ifelse(p.dat$iso3=="VEN",0,p.dat$igme_pncrisis), p.dat$year, FUN=sum)



p.dat$pnd<-rowSums(p.dat[,paste0("post",c(2:3,5:13,15:18))])




codes<- read.dta13("other_inputs/codelist region.dta")

summary(p.dat[,paste0("post",c(2:3, 5:13, 15:18))])
n.dat$whoname<- NULL
#u5<- merge(n.dat, p.dat, by=c("iso3","year"))
u5<- merge(n.dat, p.dat[,c('iso3',"year","whoreg6",paste0("post",c(2:3, 5:13, 15:18)),"pnd","lb","measout","pnmeth",
                           "malwho","mal","igme_pnd","igme_pnd_pcrisis","tb_resp","tb_nonresp","tb_resp_corr",
                           "igme_u5mr_pcrisis",
                           "tb_nonresp_corr", "pfpr","tb_resp_plus","tb_resp_plus_corr")],
           by=c("iso3", "year"))
u5<- merge(u5, codes[,c("iso3","whoname")], by="iso3")
           dim(u5)
           table(u5$year, useNA="always")
           summary(u5)
table(u5$whoreg6, useNA="always")
table(u5$whoname, useNA="always")

u5$u5d<- u5$nnd+u5$pnd
names(u5)
tapply(u5[,c("nnd")], u5$year, FUN=sum)
tapply(u5[,c("pnd")], u5$year, FUN=sum)
tapply(u5[,c("u5d")], u5$year, FUN=sum)
tapply(u5[,c("measout")], u5$year, FUN=max)
length(unique(n.dat$iso3))
length(unique(p.dat$iso3))
length(unique(u5$iso3))
#save(n.dat, p.dat, u5, file="child_cod_2000-2019_TB.RData")

names(u5)
##########################################################
#add regional and global estimates
var<- c(paste0("neo",c(2:3,5:12,19,15,17,18)),
        paste0("post",c(2:3,5:13,15:18)),
        "pnd","nnd","u5d","lb")



#####################################################################
# Estimate World
w.est<- aggregate(u5[, var], by=list(u5$year), FUN=sum)
names(w.est)[1]<- "year"
w.est$level<- "global"
#####################################################################
# Estimate Regional
reg.est<- aggregate(u5[, var], by=list(u5$whoreg6,u5$year), FUN=sum)
names(reg.est)[1:2]<- c("whoreg6","year")
reg.est$level<- "region"
u5$level<- "country"
for(v in names(u5)) {
  if(! v %in% names(w.est)) w.est[,v]<-NA
  if(! v %in% names(reg.est)) reg.est[,v]<-NA
}
temp<- rbind(w.est, reg.est, u5)
dim(u5)
dim(temp)
for(v in c(2:3,5:12,15,17,18))  {
      temp[,paste0("ufive",v)]<- temp[,paste0("neo",v)] +temp[,paste0("post",v)]
}
temp[,paste0("ufive",19)]<- temp[,paste0("neo",19)] +temp[,paste0("post",13)] +temp[,paste0("post",16)]


temp$nmr<- temp$nnd/temp$lb
temp$pnmr<- temp$pnd/temp$lb
temp$u5mr<- temp$u5d/temp$lb
temp<- temp[,c("iso3","year","level","whoreg6","whoname","nnd","pnd","u5d",
               "nmr","pnmr","u5mr","lb",
               paste0("neo",c(2:3,5:12,15,17,18,19)),
               paste0("post",c(2:3,5:13,15:18)),
               paste0("ufive",c(2:3,5:12,15,17,18,19)) ,
               "pnmeth","neometh")]



write.dta(temp, file="child_cod_2000-2019.dta")

temp<-rowSums(u5[,paste0("post",c(2:3,5:13,15:18))])
summary(temp - u5$pnd)

temp<-rowSums(u5[,paste0("neo",c(2:3,5:12,19,15,17:18))])
summary(temp - u5$nnd)

est<- u5
summary(est)
