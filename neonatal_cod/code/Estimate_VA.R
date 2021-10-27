
#### Analysis of samples
library(R2jags)
library(doParallel)
library(coda)
library(lattice)
library(dplyr)


load(file="../data/data_missing_2020_05_27.RData") ## Load data
source("functions_VA.R")        ## load functions



### check distribution of variables for the model:
histogram(studies$u5mr)
histogram(studies$nmr)
histogram(studies$dpt)
histogram(studies$gfr)
histogram(studies$lbwrate)
histogram(studies$sba)
histogram(studies$bcg)
histogram(studies$pab)
histogram(studies$femlit)

# For sba make sure is all between 0 and 1
studies <- mutate(studies, sba=if_else(sba>1, sba/100, sba))

## Values for parameters in the functions
my.modl  <- "LASSO_final_VA.txt"    # File with model
my.lamb  <- 10
my.rsdl  <- 0.07
my.quad  <- FALSE
my.iter  <- 10000  # number of mcmc SAMPLES
my.burn  <- 5000   # number of draws used for burn-in
my.chas  <- 4      # number of chains
my.thin  <- 1      # thin factor
my.boot  <- 100    # number of times bootstrap samples
my.prin  <- 1      # don't print results while calculating
my.samc  <- T      # Save JAGS object
my.sade  <- F
my.savx  <- F
my.vxr   <- vxr
my.vxf   <- vxf
my.vdt   <- vdt

lams <- c(225, 125)
resd <- c(0.14, 0.21)

####################################
#######   MODEL ESTIMATION   #######
####################################

for (i in 1:length(lams)){
  nam <- paste("va",lams[i],resd[i],sep="_")
  assign(nam,f.e1(STUD=studies, DEAT=deaths, LAMB=lams[i], RSDL=resd[i], SAVX=F, SADE=F, NAME=paste("va",lams[i],resd[i],sep="_")))
  save(list=paste("va",lams[i],resd[i],sep="_"), file = paste0("../results/VA_",lams[i],"_",resd[i],".RData"))
  do.call(rm, list(paste("va",lams[i],resd[i],sep="_")))
}




