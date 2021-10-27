
#############################
# admin for Amy
#setwd("H:/MCEE/analysis")
rm(list=ls())
Start.time<- Sys.time()

#############################
#### Analysis of samples
library(sirt)
library(R2jags)
library(doParallel)
library(coda)
library(lattice)



## Load data
load(file="../data/PN_VR_data.RData")
vxf #covariates
vxr #Random effect identifier
vdt #causes
dim(studies)





#Written by Amy Mulick and David Prieto
source("functions_VRMod.R")


test.label<- "Test0"


vxf<-  c( "Iwhoreg6_3",
         "year",
         "ln_5q0_cov",
         "ln_ppp_cov",
         "water_cov","Hib3_cov","skilled_ba_cov")




############  Esimate Models with apparent precission (all data) ################
for(lam in c(5)) {
for(sd in c(0.07)) {
  vxf



  Start<- Sys.time()
  print(Start)
  cat("\n\n LAMBDA ",lam," , sd = ", sd, "  \n\n")


  out<- f.e1(PARA=1,        #1 - use parallel computing, 0 - use serial computing
             STUD=studies,  #Input study data
             DEAT=deaths,   #Deaths from input studies (in misclassification matrix)
             VDT=vdt,       #Causes
             VXF=vxf,       #Covariates
             VXR=vxr,       #Random Effect identifier
             MODL="LASSO_VR.txt", #Text file with JAGS model defintion
             LAMB=lam,      #Lambda value for Lasso
             RSDL=sd,       #For prior of standard deviation of random effects, maximum of SD
             QUAD=0,        #NULL if no quadratic terms desired
             CHAS=4,        #Number of MCMC chains
             ITER=30000,      #Number of MCMC iterations
             BURN=23000,        #Number of iterations for burn-in of MCMC
             THIN=5,        #Thining factor for MCMC
             SEED=77745123,      #random seed for jags.parallel (for jags use set.seed() first)
             PRIN=TRUE,     #Print computation stage
             SAMC=TRUE,     #T/F save MCMC samples?
             SAVX=TRUE,     #Save explanatory variables
             SADE=TRUE,     #Save death matrix
             NAME=NULL      #Optional model names
             )


  #names(out)
  #out$param
  #names(out$output$BUGSoutput)

  print(Sys.time() - Start)

  #names(jm.r0)
  #head(summary(jm.r0))

  ######################################################
  #convergence statistics
  class(out$output$BUGSoutput)
  #summ <-  sirt::mcmc.list.descriptives( as.mcmc(out$output$BUGSoutput) )
  #warnings()
  summ<-as.data.frame(out$output$BUGSoutput$summary)
    names(summ)
    cat("\n\nSummary Rhat\n")
  print(summary(summ$Rhat))
  #hist(summ$Rhat, main=paste("Lambda resampled, sd max = ",sd))

  eb0   <- out$output$BUGSoutput$median$B # fixed effects
  cat("\n\n BETA estimate\n\n")
  print(eb0)

  save(out,
       #summ, jm.r0,
       #lambda, sd,
       #studies, data.predict, deaths,
       #vxf, #predictors
       #vdt, #causes
       file=paste0("../results/VR.RData"))

} #end loop over sds
} #end loop over lambdas



#summ
Sys.time() - Start.time
