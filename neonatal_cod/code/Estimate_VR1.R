#############################
# admin
rm(list=ls())
#############################

#### Analysis of samples
library(R2jags)
library(doParallel)
library(coda)
library(lattice)


## Values for parameters in the functions
my.func  <- "Analysis_functions_CV_9.R" # File with functions
my.data  <- "../data/data1_VR_2020_05_18.RData"          # File with R data
my.modl  <- "LASSO_final_VR.txt"    # File with model
my.lamb  <- 10
my.rsdl  <- 0.07
my.quad  <- FALSE
my.iter  <- 1000  # number of mcmc SAMPLES
my.burn  <- 500   # number of draws used for burn-in
my.chas  <- 4      # number of chains
my.thin  <- 1      # thin factor
my.prin  <- 1      # don't print results while calculating
my.samc  <- T      # Save JAGS object
my.sade  <- F
my.savx  <- F

load(file=paste(my.data))    ## Load data
source(paste(my.func))       ## load functions

my.vxr   <- vxr
my.vxf   <- vxf
my.vdt   <- vdt

lams <- c(325)
l.lams <- length(lams)
resd <- seq(0.07, 0.07, 0.07)
l.resd <- length(resd)



####################################
#######   MODEL ESTIMATION   #######
####################################

#setwd("~/Dropbox/Amy Parallel/VR_19.05.2020/model_estimation")
for (i in lams){
  nam <- paste("vr1",i,sep="_")
  assign(nam,f.e1(STUD=studies1, DEAT=deaths1, LAMB=i, SAVX=F, SADE=F, NAME=paste("vr1",i,sep="_")))
  save(list=paste("vr1",i,sep="_"), file = paste0("../results/VR1.RData"))
  do.call(rm, list(paste("vr1",i, sep="_")))
}



####################################
#######   MODEL CHECKING   ########
####################################

#for (i in lams) load(paste0("../results/VR1.RData"))
#vr1nams <- ls()[grep("vr1_", ls())]
#mod.names <- rownames(get(vr1nams[1])[["output"]][["BUGSoutput"]][["summary"]])

#for (v in vr1nams) {
#  pdf(file = paste0("model_check_",v,".pdf"), paper="a4")
#  plot(as.mcmc(get(v)$output)[][,grep("B", mod.names)],density = F, sub = paste0("model==",v))
#  gelman.plot(as.mcmc(get(v)$output)[][,grep("B", mod.names)], ylim=c(1,2))
#  dev.off()
#}


