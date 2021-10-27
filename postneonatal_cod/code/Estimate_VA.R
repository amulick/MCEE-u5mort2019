
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
load(file="../data/PN_VA_data.RData")
vxf #covariates
vxr #Random effect identifier
vdt #causes

names(studies)
source("functions_VA.R")


vxf <- c("Iagelower_1","u5mr","pfpr","meningitis_epi","gni","mcv","yr_mid","underwt","sanitation")
vxf



############  Esimate Models with apparent precission (all data) ################
for(lam in c(10)) {
for(sd in c(0.14)) { #c(0.21, 0.28, 0.35, 0.42 )) {


  Start<- Sys.time()
  print(Start)
  cat("\n\n LAMBDA ",lam," , sd = ", sd, "  \n\n")


  out<- f.e1(PARA=1,        #1 - use parallel computing, 0 - use serial computing
             STUD=studies,  #Input study data
             DEAT=deaths,   #Deaths from input studies (in misclassification matrix)
             VDT=vdt,       #Causes
             VXF=vxf,       #Covariates
             VXR=vxr,       #Random Effect identifier
             MODL="LASSO_VA.txt", #Text file with JAGS model defintion
             LAMB=lam,      #Lambda value for Lasso
             RSDL=sd,       #For prior of standard deviation of random effects, maximum of SD
             QUAD=0,        #NULL if no quadratic terms desired
             CHAS=4,        #Number of MCMC chains
             ITER=15000,      #Number of MCMC iterations
             BURN=5000,        #Number of iterations for burn-in of MCMC
             THIN=5,        #Thining factor for MCMC
             SEED=123,      #random seed for jags.parallel (for jags use set.seed() first)
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
       file=paste0("../results/VA.RData"))

} #end loop over sds
} #end loop over lambdas




Sys.time() - Start.time
