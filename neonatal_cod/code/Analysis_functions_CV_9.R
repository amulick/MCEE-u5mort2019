
## Functions needed for Bootstrap and analysis



##################################################Z
####
####   Model estimation  (one model)
####  


f.e1 <- function(PARA=1, STUD=my.stud, DEAT=my.deat, VDT=my.vdt, 
                 VXF=my.vxf, VXR=my.vxr, MODL=my.modl, LAMB=my.lamb, 
                 RSDL=my.rsdl, QUAD=my.quad, CHAS=my.chas, ITER=my.iter, 
                 BURN=my.burn, THIN=my.thin, SEED=123, PRIN=my.prin, 
                 SAMC=my.samc, SAVX=my.savx, SADE=my.sade, NAME=NULL) {
  ##
  ##  Function for ENTIRE sample
  ##
  ##  PARA = 1 to run CHAINS in parallel and 0 to run CHAINS in serial
  ##  ------------------------- #
  ##  STUD = dataset of studies
  ##  DEAT = dataset of reported causes of death by study with missclassification matrix
  ##  VDT  = names of true causes of death (string vector)
  ##  VXF  = names of fixed effects variables (string vector)
  ##  VXR  = name of random effects variable (string vector)
  ##  -------------- #
  ##  MODL = txt file with BUGS model
  ##  LAMB = Lambda value for Poisson
  ##  RSDL = limit for the random effects SD
  ##  QUAD = NULL if no quadratic terms desired (default); !is.null(QUAD) adds quadratic terms
  ##  -------------- #
  ##  PARA = 0 run CHAINS in serial, PARA!=0 run CHAINS in parallel
  ##  CHAS = number of chains
  ##  ITER = number of mcmc SAMPLES
  ##  BURN = number of draws used for burn-in (currently arbitrarily set to 500)
  ##  THIN = thin factor 
  ##  SEED = random seed for jags.parallel (for jags use set.seed() first)
  ##  -------------- #
  ##  PRIN = Print stage of computation
  ##  SAMC = Save MCMC samples? TRUE/FALSE
  ##  SAVX = T, Save explanatory variables
  ##  SADE = T, Save matrix of deaths
  ##  NAME = Optional name for the model
  ##
  ptime <- proc.time()
  ## Prepare DATA
  if (QUAD) {
    STUD[, paste0(VXF, ".q")] <- STUD[, VXF]^2
    VXF <- c(VXF, paste0(VXF, ".q"))
  }
  # Add numeric variable for Re group to data
  STUD$rG <- as.numeric(as.factor(STUD[,VXR]))
  # drop deaths from studies not in estimation sample
  DEAT <- droplevels(DEAT[DEAT$sid %in% STUD$sid, ])
  # Save DEAT if required for return
  if(SADE) sde <- DEAT else sde <- NULL
  
  ## Prepare JAGS DATA
  # Constants for model
  K  <- length(VXF)+1            # Number of fixed-effects parameter
  S  <- dim(STUD)[1]             # Number of studies 
  C  <- length(VDT)              # Number of true causes of DEAT
  rN <- length(unique(STUD$rG))  # number of RE groups
  # First reported death
  fi <- as.numeric(tapply(1:nrow(DEAT), DEAT$sid, min))
  # Last reported death
  la <- as.numeric(tapply(1:nrow(DEAT), DEAT$sid, max))
  # Number of deaths in the study
  N <- as.numeric(tapply(DEAT$n, DEAT$sid, sum))
  # Create GM matrix
  GM <- as.matrix(DEAT[, c(VDT,"n")])
  ## Scale & centre continuos covariates (including squared terms)
  XM <- cbind(CONS=rep(1,S), as.matrix(STUD[,VXF]))
  XM[, -1] <- apply(XM[, -1], 2, scale)
  # List with JAGS data
  jd.r1 <- list(S = S, C = C, K = K, rN = rN, lambda = LAMB, 
                rsdlim = RSDL, fi = fi, la = la, N = N,
                rG = STUD$rG, XM = XM, GM = GM)
  
  ## Parameters to follow
  jp.r1 <- c("B","re","sd")
  ## Initial values, randomly generated from distributions with large variance
  ## except 'sd' parameters, which have fixed starting values
  ji.r1 <- function() list(B=matrix(c(rep(NA,K),rnorm((K*(C-1)),0,8)),K,C),
                           re=matrix(c(rep(NA,rN),rnorm((rN*(C-1)),0,8)),rN,C),
                           sd=c(NA,rep(0.03,(C-1))))
  
  # Remove all unnecessary objects before calling JAGS
  rm(DEAT, XM, GM, la, fi, N)
  # END NEW
  
  ## SAMPLES
  if(PARA==1){
    ## Set parameters in the Global environment
    assign("my.iter", ITER, envir=.GlobalEnv)
    assign("my.burn", BURN, envir=.GlobalEnv)
    assign("my.thin", THIN, envir=.GlobalEnv)
    jm.r1 <- jags.parallel(data=jd.r1, parameters.to.save=jp.r1, 
                           inits=ji.r1, model.file=paste(MODL), 
                           jags.seed = SEED, n.chains=CHAS, 
                           n.iter=my.iter, n.burnin=my.burn, 
                           n.thin=my.thin, 
                           export_obj_names=c("my.iter","my.burn","my.thin")) # Parallel chains
  }
  if(PARA==0){
    jm.r1 <- jags(data=jd.r1, parameters.to.save=jp.r1, inits=ji.r1, 
                  model.file=paste(MODL), n.chains=CHAS, n.iter=ITER,
                  n.burnin=BURN, n.thin=THIN) # NON-parallel chains
  }
  
  ###  Prepare Return list  ###
  # Define parameters of the model
  para <- list(name = NAME, time = ptime, PARA = PARA, 
               VDT = VDT, VXF = VXF, VXR = VXR, MODL = MODL, 
               LAMB = LAMB, RSDL = RSDL, QUAD = QUAD)
  # Save means and SD of predictors for standardisation:
  dv <- data.frame(xvar=VXF, mean = colMeans(STUD[, VXF], na.rm = T), 
                   sd = apply(STUD[, VXF], 2, sd, na.rm=T), 
                   stringsAsFactors = F, row.names = NULL)
  # decide what prediction variables to save
  if(SAVX==T) svx <- c("sid",VXR,"rG",VXF) else svx <- c("sid",VXR,"rG")
  # Remove MCMC samples if not wanted
  if(SAMC==F){
    jm.r1$BUGSoutput$sims.array <- NULL
    jm.r1$BUGSoutput$sims.list <- NULL # don't save these
    jm.r1$BUGSoutput$sims.matrix <- NULL # don't save these
  }
  ptime <- proc.time() - ptime
  if(PRIN==1) print(paste("time = ",ptime[3]))
  return(list(param = para, Vars = dv, Studies = STUD[,svx], 
              Deaths = sde, output=jm.r1))
}

###
###  END OF FUNCTION
###
################################################ #




################################################ #
####
####   Prediction from ONE model estimation in other dataset
####  

f.p1 <- function(MO, STUD, DEAT=NULL, ADPA=NULL, PFIX=T, PRAN=T, REME=0, PSUF=NULL){
  ##
  ##  MO   = model to use
  ##  STUD = data from studies to predict
  ##  DEAT = deaths in studies to predict
  ##  ADPA = Add estimation parameters to the table
  ##  PFIX = Predict with fixed effects?
  ##  PRAN = Predict with fixed+Random effects?
  ##  REME = put "0" in Re for studies without a Re in from estimation (otherwise use observed mean of Re)
  ##  PSUF = suffix for prediction  
  
  S  <- dim(STUD)[1]        # number of studies in estimation set
  K  <- dim(MO$Vars)[1]     # number of variables in model
  TC <- MO$param$VDT        # vector of true causes of death

  
  # prepare raw variables dataset from prediction sample
  DX <- cbind(rep(1,S), as.matrix(STUD[,MO$Vars$xvar]))
  meanMat <- matrix(MO$Vars$mean, nrow = nrow(DX), ncol = K, byrow = T)
  sdMat <- matrix(MO$Vars$sd, nrow = nrow(DX), ncol = K, byrow = T)
  DX[,-1] <- (DX[,-1] - meanMat)/sdMat
  rm(meanMat, sdMat)

  # PREDICTION True causes distribution WITH FIXED EFFECTS
  if(PFIX==T){
    PF <- DX %*% MO$output$BUGSoutput$mean$B    # linear predictor
    PF <- cbind(rep(1,S), exp(PF))  # OdSTUD 
    PF <- cbind(sid=STUD[,c("sid","lab")], as.data.frame(PF/rowSums(PF))) # Proportions
    names(PF) <- c("sid","lab",TC)             # Put names on variables
  }

  # PREDICT True causes distribution WITH RANDOM EFFECTS
  if(PRAN==T){
    STUD <- merge(STUD, unique(MO$Studies[,c(2,3)]), by=names(MO$Studies)[2], all.x=T, all.y=F)
    # Prepare average andom effects
    PR <- MO$output$BUGSoutput$mean$re[STUD$rG, ]   # Prepare RE matrix
    MR <- colMeans(MO$output$BUGSoutput$mean$re, na.rm = T)    # unweighted means of RE
    if(REME == 0) MR <- rep(0, length(MR))
    PR[is.na(PR)] <- MR[ceiling(which(is.na(PR))/nrow(PR))] # Put MR if no RE for study
    PR <- DX %*% MO$output$BUGSoutput$mean$B + PR # linear PRedictor
    PR <- cbind(rep(1,S), exp(PR))   # OdSTUD 
    PR <- cbind(sid=STUD[,c("sid","lab")], as.data.frame(PR/rowSums(PR))) # PRoportions
    names(PR) <- c("sid","lab",TC)   # Put names on variables
  }
  
  ## PREDICT OBSERVED DEATHS (If data available)
  if(!is.null(DEAT)){
        if (PFIX | PRAN) {
      # Count number of appearences of each sid
      nsid <- as.vector(table(DEAT$sid))
      # Vector with number of cases in each sid
      nn <- as.vector(tapply(DEAT$n, DEAT$sid, sum))
      nn <- rep(nn, nsid)
      # fixed effects predictions
      if (PFIX){
        PFmat <- PF[rep(seq_len(nrow(PF)), times = nsid), ]
        PFmat <- PFmat[order(PFmat$sid), ]
        DEAT[, paste0("pf",PSUF)] <- nn * apply(DEAT[, TC] * PFmat[, TC], 1, sum)
      }
      # random effects predictions
      if (PRAN){
        PRmat <- PR[rep(seq_len(nrow(PR)), times = nsid), ]
        PRmat <- PRmat[order(PRmat$sid), ]
        DEAT[, paste0("pr",PSUF)] <- nn * apply(DEAT[, TC] * PRmat[, TC], 1, sum)
      }
    }
    # Was study in estimation sample?
    DEAT[, paste0("e.samp",PSUF)] <- DEAT$sid %in% MO$Studies$sid
    DEAT$p.name <- MO$param$name
    for(v in ADPA) DEAT[,paste0("p.",v)] <- MO$param[[v]]
  }
  rtn <- list()
  rtn[["Para_pre"]] <- match.call()
  rtn[["Para_est"]] <- MO$param
  if(!is.null(PF)) rtn[["Pfix"]] <- PF
  if(PRAN==T)      rtn[["Pran"]] <- PR
  if(!is.null(DEAT)) rtn[["Pdat"]] <- DEAT
  return(rtn)
}

###
###  END OF FUNCTION
###
#################################################Z



################################################Z
####
####   k-fold ESTIMATIONS 
####  


f.cv <- function(CVNA, KFOL, STUD, DEAT, VDT=my.vdt, VXF=my.vxf, 
                 VXR=my.vxr, LAMB=my.lamb, RSDL=my.rsdl, QUAD=my.quad, 
                 MODL=my.modl, CHAS=my.chas, ITER=my.iter, BURN=my.burn,
                 THIN=my.thin, PRIN=my.prin, SAMC=my.samc, SAVX=my.savx, 
                 SADE=my.sade, CORE=KFOL, JS=123, PS=98765){
  ##  CVNA = Name of cross validation set
  ##  KFOL = Number of CV partitions of the data
  ##  ---------------------------------------- #
  ##  STUD = dataset of studies
  ##  DEAT = dataset of reported causes of death by study with missclassification matrix
  ##  VDT  = names of true causes of death (string vector)
  ##  VXF  = names of fixed effects variables (string vector)
  ##  VXR  = name of random effects variable (string vector)
  ##  -------------- #
  ##  LAMB = Lambda value for Poisson
  ##  RSDL = limit for the random effects SD
  ##  QUAD = NULL if no quadratic terms desired (default); !is.null(QUAD) adds quadratic terms
  ##  -------------- #
  ##  MODL = txt file with BUGS model
  ##  PARA = 0 run CHAINS in serial, PARA!=0 run CHAINS in parallel
  ##  CHAS = number of mcmc SAMPLES
  ##  ITER = number of draws used for burn-in (currently arbitrarily set to 500)
  ##  BURN = number of chains
  ##  THIN = thin factor 
  ##  --------------- #
  ##  PRIN = Print stage of computation
  ##  SAMC = Save JAGS object for MCMC diagnostics? TRUE/FALSE
  ##  SAVX = T, Save explanatory variables
  ##  SADE = T, Save matrix of deaths
  ##  CORE = Number of cores to use for CV samples (>1 = parallelise CV partitions (not chains in mcmc (default)), 1=parallel chains (not CV partitions), 0= do not parallelise anything)
  ##  --------------- #
  ##  JS   = JAGS seed, i.e. the random seed for MCMC sampling, default 123
  ##  PS   = Partition seed, i.e. the random seed for generating k-fold partitions, default 98765
  
  ## Set parameters in the Global environment
  assign("my.deat", DEAT, envir=.GlobalEnv)  
  assign("my.vdt",  VDT,  envir=.GlobalEnv) 
  assign("my.vxf",  VXF,  envir=.GlobalEnv) 
  assign("my.vxr",  VXR,  envir=.GlobalEnv) 
  assign("my.modl", MODL, envir=.GlobalEnv)
  assign("my.lamb", LAMB, envir=.GlobalEnv)
  assign("my.rsdl", RSDL, envir=.GlobalEnv)
  assign("my.quad", QUAD, envir=.GlobalEnv)
  assign("my.chas", CHAS, envir=.GlobalEnv)
  assign("my.iter", ITER, envir=.GlobalEnv)
  assign("my.burn", BURN, envir=.GlobalEnv)
  assign("my.thin", THIN, envir=.GlobalEnv)
  assign("my.prin", PRIN, envir=.GlobalEnv)
  assign("my.samc", SAMC, envir=.GlobalEnv)
  assign("my.savx", SAVX, envir=.GlobalEnv)
  assign("my.sade", SADE, envir=.GlobalEnv)
  assign("my.cvna", CVNA, envir=.GlobalEnv)
  
  ## Create the K-fold partition
  # STRATIFIED SAMPLING BY RANDOM EFFECT TERM
  S1 <- length(unique(STUD[,VXR]))               # number of strata for sampling 
  set.seed(PS)
  GCV  <- sample(rep(c(1:KFOL), S1)[1:S1], S1)   # create a vector of cross validation groups, randomly ordered
  names(GCV) <- as.character(unique(STUD[,VXR])) # name vector according to (ordered) names of random effects
  STUD$GCV  <- GCV[as.character(STUD[,VXR])]     # add cross validation group to dataset
  
  assign("my.stud", STUD, envir=.GlobalEnv)  # Put study data in glopbal environnment
  
  ptime <- proc.time()
  if(CORE<2){
    STA <- list()
    for(i in 1:KFOL) STA[[k]] <- f.e1(PARA = CORE, STUD=STUD[STUD$GCV!=i,],
                                      NAME = paste(CVNA,i,sep="."), SEED = JS)
  } else {  
    nc <- detectCores()
    cl <- makeCluster(min(CORE, nc))
    registerDoParallel(cl)
    STA <- foreach(i = c(1:KFOL), .packages = c("R2jags"),
                   .export=c("f.e1", "my.deat", "my.vdt", "my.vxf", "my.vxr", "my.modl",
                             "my.lamb", "my.rsdl", "my.quad", "my.chas", "my.iter",
                             "my.burn", "my.thin", "my.prin", "my.samc", "my.savx",
                             "my.sade")) %dopar% f.e1(PARA = 0,
                                                      STUD = STUD[STUD$GCV != i, ],
                                                      NAME = paste(CVNA, i, sep="."))
    stopCluster(cl)
  }
  
  ptime <- proc.time() - ptime
  print(paste("Time of total CV estimation:", ptime[3]))
  names(STA) <- paste(CVNA,c(1:KFOL),sep=".")
  assign(CVNA, STA , envir=.GlobalEnv)
}


###
###  END OF FUNCTION
###
################################################ #




################################################ #
####
####   Prediction from a set of CV results
####  

f.pcv <- function(CVS, STUD, DEAT, ADPA=NULL, PFIX=T, PRAN=T, REME=0){
  ##
  ##  CVs  = list of text names of models to use
  ##  STUDs = data from studies to predict
  ##  DEATs = deaths in studies to predict
  ##  ADPAs = Add estimation parameters to the table
  ##  PFIXs = Predict with fixed effects?
  ##  PRANs = Predict with fixed+Random effects?
  ##  REMEs = put "0" in Re for studies without a Re in from estimation (otherwise use observed mean of Re)
  
  for(cv in CVS){
    nc <- detectCores()
    cl <- makeCluster(min(length(cv), nc))
    registerDoParallel(cl)
    pcv <- foreach(M=iter(get(cv)), .combine=rbind, .export=c("f.p1")) %dopar% f.p1(M, STUD, DEAT, ADPA, PFIX, PRAN, REME)$Pdat
    stopCluster(cl)
    if(cv==CVS[1]) pcvs <- pcv else pcvs <- rbind(pcvs, pcv)
  }
  return(pcvs)
}

###
###  END OF FUNCTION
###
#################################################Z






##################################################Z
####
####   Summarysing a table of predictions
####  

f.sp <- function(PT, FUNC, vt=NULL){
  ##
  ##  PT  = prediction table
  ##  vt  = variables to stratify summaries
  PT$All <- "All"
  if(is.null(vt)) vt <- "All"
  lvt <- list()
  for(i in 1:length(vt)) lvt[[i]] <- PT[,vt[i]]
  tt <- aggregate(PT[,c("n")], lvt, sum) # count number of deaths in strata
  names(tt) <- c(vt, "N")

  if(is.null(PT$pf)==F){
    PT$pf1 <- (PT$pf - PT$n)^2
    PT$pf2 <- (PT$pf1)*PT$pf
    trf <- aggregate(PT[,c("pf2")], lvt, sum)
    names(trf) <- c(vt,"pf2")
    tt <- merge(tt,trf, by=vt)
    tt$MSE.f <- tt$pf2/tt$N
    tt$RMSE.f <- sqrt(tt$pf2/tt$N)
    for(c in c("pf2", "MSE.f", "RMSE.f")) tt[,c] <- formatC(tt[,c], 1, format="f")
  }

  if(is.null(PT$pr)==F){
    PT$pr1 <- (PT$pr - PT$n)^2
    PT$pr2 <- (PT$pr1)*PT$pr
    trr <- aggregate(PT[,c("pr2")], lvt, sum)
    names(trr) <- c(vt,"pr2")
    tt <- merge(tt,trr, by=vt)
    tt$MSE.r <- tt$pr2/tt$N
    tt$RMSE.r <- sqrt(tt$pr2/tt$N)
    for(c in c("pr2", "MSE.r", "RMSE.r")) tt[,c] <- formatC(tt[,c], 1, format="f")
  }
  return(tt)
}

###
###  END OF FUNCTION
###
#################################################Z






##################################################Z
####
####   Model estimation  SAVING CHAINS 
####  


f.es <- function(PARA=1,       STUD=my.stud, DEAT=my.deat, VDT=my.vdt, 
                 VXF=my.vxf,   VXR=my.vxr,   MODL=my.modl, LAMB=my.lamb, 
                 RSDL=my.rsdl, QUAD=my.quad, CHAS=my.chas, BURN=my.burn, 
                 THIN=my.thin, SAIT=my.sait, SAFI=my.safi, SEED=123, PRIN=my.prin, 
                 SAMC=my.samc, SAVX=my.savx, SADE=my.sade, NAME) {
  ##
  ##  Function for ENTIRE sample
  ##
  ##  PARA = way to run and files to save: 0=serial, 1=parallel, >1 = parallel+save in PARA files
  ##  ------------------------- #
  ##  STUD = dataset of studies
  ##  DEAT = dataset of reported causes of death by study with missclassification matrix
  ##  VDT  = names of true causes of death (string vector)
  ##  VXF  = names of fixed effects variables (string vector)
  ##  VXR  = name of random effects variable (string vector)
  ##  -------------- #
  ##  MODL = txt file with BUGS model
  ##  LAMB = Lambda value for Poisson
  ##  RSDL = limit for the random effects SD
  ##  QUAD = NULL if no quadratic terms desired (default); !is.null(QUAD) adds quadratic terms
  ##  -------------- #
  ##  PARA = 0 run CHAINS in serial, PARA!=0 run CHAINS in parallel
  ##  CHAS = number of mcmc chains
  ##  BURN = number of draws used for burn-in (currently arbitrarily set to 500)
  ##  THIN = thin factor 
  ##  SAIT = number of draws to finally save per chain per file
  ##  SEED = random seed for jags.parallel (for jags use set.seed() first)
  ##  -------------- #
  ##  PRIN = Print stage of computation
  ##  SAMC = Save MCMC samples? TRUE/FALSE
  ##  SAVX = T, Save explanatory variables
  ##  SADE = T, Save matrix of deaths
  ##  NAME = COMPULSORY name for the model
  ##
  
  ## Prepare DATA
  if (QUAD) {
    STUD[, paste0(VXF, ".q")] <- STUD[, VXF]^2
    VXF <- c(VXF, paste0(VXF, ".q"))
  }
  # Add numeric variable for Re group to data
  STUD$rG <- as.numeric(as.factor(STUD[,VXR]))
  # drop deaths from studies not in estimation sample
  DEAT <- droplevels(DEAT[DEAT$sid %in% STUD$sid, ])
  # Save DEAT if required for return
  if(SADE) sde <- DEAT else sde <- NULL
  
  ## Prepare JAGS DATA
  # Constants for model
  K  <- length(VXF)+1            # Number of fixed-effects parameter
  S  <- dim(STUD)[1]             # Number of studies 
  C  <- length(VDT)              # Number of true causes of DEAT
  rN <- length(unique(STUD$rG))  # number of RE groups
  ITER = THIN*SAIT   # number of iterations before thin
  # First reported death
  fi <- as.numeric(tapply(1:nrow(DEAT), DEAT$sid, min))
  # Last reported death
  la <- as.numeric(tapply(1:nrow(DEAT), DEAT$sid, max))
  # Number of deaths in the study
  N <- as.numeric(tapply(DEAT$n, DEAT$sid, sum))
  # Create GM matrix
  GM <- as.matrix(DEAT[, c(VDT,"n")])
  ## Scale & centre continuos covariates (including squared terms)
  XM <- cbind(CONS=rep(1,S), as.matrix(STUD[,VXF]))
  XM[, -1] <- apply(XM[, -1], 2, scale)
  # List with JAGS data
  jd.r1 <- list(S = S, C = C, K = K, rN = rN, lambda = LAMB, 
                rsdlim = RSDL, fi = fi, la = la, N = N,
                rG = STUD$rG, XM = XM, GM = GM)
  
  ## Parameters to follow
  jp.r1 <- c("B","re","sd")
  ## Initial values, randomly generated from distributions with large variance
  ## except 'sd' parameters, which have fixed starting values
  ji.r1 <- function() list(B=matrix(c(rep(NA,K),rnorm((K*(C-1)),0,8)),K,C),
                           re=matrix(c(rep(NA,rN),rnorm((rN*(C-1)),0,8)),rN,C),
                           sd=c(NA,rep(0.03,(C-1))))
  
  # Remove all unnecessary objects before calling JAGS
  rm(DEAT, XM, GM, la, fi, N)
  # END NEW
  
  ## SAMPLES
  ptime <- proc.time()
  if(PARA==0){  # SEQUENTIAL run
    jm.r1 <- jags(data=jd.r1, parameters.to.save=jp.r1, inits=ji.r1, 
                  model.file=paste(MODL), n.chains=CHAS, n.iter=ITER,
                  n.burnin=BURN, n.thin=THIN) # NON-parallel chains
  }
  if(PARA==-1){  # PARALLEL
    ## Set parameters in the Global environment
    assign("my.iter", ITER, envir=.GlobalEnv)
    assign("my.burn", BURN, envir=.GlobalEnv)
    assign("my.thin", THIN, envir=.GlobalEnv)
    jm.r1 <- jags.parallel(data=jd.r1, parameters.to.save=jp.r1, 
                           inits=ji.r1, model.file=paste(MODL), 
                           jags.seed = SEED, n.chains=CHAS, 
                           n.iter=my.iter, n.burnin=my.burn, 
                           n.thin=my.thin, 
                           export_obj_names=c("my.iter","my.burn","my.thin")) # Parallel chains
  }
  if(PARA>0){   # PARALLEL + saveJAGS
    jm.r1 <- saveJAGS(data=jd.r1, inits=ji.r1, params=jp.r1, modelFile=paste(MODL),
                      chains=CHAS, sample2save=SAIT, nSaves=PARA, burnin=BURN, thin=THIN,
                      fileStub=NAME )
  }
  ptime <- proc.time() - ptime
  if(PRIN==1) print(paste("time = ",ptime[3]))
  ###  Prepare Return list  ###
  # Define parameters of the model
  para <- list(name = NAME, time = ptime, PARA = PARA, 
               VDT = VDT, VXF = VXF, VXR = VXR, MODL = MODL, 
               LAMB = LAMB, RSDL = RSDL, QUAD = QUAD)
  # Save means and SD of predictors for standardisation:
  dv <- data.frame(xvar=VXF, mean = colMeans(STUD[, VXF], na.rm = T), 
                   sd = apply(STUD[, VXF], 2, sd, na.rm=T), 
                   stringsAsFactors = F, row.names = NULL)
  # decide what prediction variables to save
  if(SAVX==T) svx <- c("sid",VXR,"rG",VXF) else svx <- c("sid",VXR,"rG")
  # Remove MCMC samples if not wanted
  if(SAMC==F){
    jm.r1$BUGSoutput$sims.array <- NULL
    jm.r1$BUGSoutput$sims.list <- NULL
    jm.r1$BUGSoutput$sims.matrix <- NULL
  }
  return(list(param = para, Vars = dv, Studies = STUD[,svx], 
              Deaths = sde, output=jm.r1))
}

###
###  END OF FUNCTION
###
################################################ #






################################################ #
####
####   BETTER PREDICTION FUNCTIONS by AMY
####  


f.p2 <- function(MO, STUD, ID=c("isocode", "year"), DEAT=NULL, ADPA=NULL, PFIX=T, PRAN=T, REME=0, PSUF=NULL, PERIOD='both'){
  ##  MO   = model to use
  ##  STUD = data from studies to predict
  ##  ID   = variables that uniquely identify observations
  ##  DEAT = deaths in studies to predict
  ##  ADPA = Add estimation parameters to the table
  ##  PFIX = Predict with fixed effects?
  ##  PRAN = Predict with fixed+Random effects?
  ##  REME = put "0" in Re for studies without a Re in from estimation (otherwise use observed mean of Re)
  ##  PSUF = suffix for prediction
  ##  PERIOD = 'early', 'late' or 'both' (only applies to VA model)
  
  S  <- dim(STUD)[1]        # number of studies in estimation set
  K  <- dim(MO$Vars)[1]     # number of variables in model
  TC <- MO$param$VDT        # vector of true causes of death
  LA <- MO$param$LAMB       # lambda in model (added 6 Nov 2019)
  RE <- MO$param$RSDL       # RE-SD in model (added 6 Nov 2019)
  
  # Add random effect group to the prediction data
  if (PRAN) {
    if (names(MO$Studies)[2] %in% names(STUD)) {
      STUD$order <- 1:dim(STUD)[1]                    #added 6 Nov 2019
      STUD <- merge(STUD, unique(MO$Studies[,c(2,3)]), 
                    by=names(MO$Studies)[2], all.x=T, all.y=F)
      groupcat <- names(MO$Studies)[2]
      STUD <- STUD[order(STUD$order), ]               #added 6 Nov 2019
    } else {
      PRAN <- F
      # If user wants PRAN, but not possible, at least PFIX
      PFIX <- T
      warning(paste0("Random effects not included because grouping category '",
                     names(MO$Studies)[2], "' is not available in the prediction dataset."))
    }
  }
  
  if (all(c('per.early', 'per.late') %in% MO$param$VXF) &
      !all(c('per.early', 'per.late') %in% names(STUD))) {
    STUD$per.early <- 0
    STUD$per.late <- 0
    if (PERIOD == 'both') {
      # Mean value of these two 'covariates' set to 0.
      MO$Vars$mean[1:2] <- 0
    }
    if (PERIOD == 'early') {
      STUD$per.early <- 1
      # Mean value for 'late' set to 0
      MO$Vars$mean[2] <- 0
    }
    if (PERIOD == 'late') {
      STUD$per.late <- 1
      # Mean value for 'early' set to 0
      MO$Vars$mean[1] <- 0
    }
  }
  
  # prepare raw variables dataset from prediction sample
  DX <- cbind(rep(1,S), as.matrix(STUD[,MO$Vars$xvar]))
  meanMat <- matrix(MO$Vars$mean, nrow = nrow(DX), ncol = K, byrow = T)
  sdMat <- matrix(MO$Vars$sd, nrow = nrow(DX), ncol = K, byrow = T)
  DX[,-1] <- (DX[,-1] - meanMat)/sdMat
  rm(meanMat, sdMat)
  
  # PREDICTION True causes distribution WITH FIXED EFFECTS
  if (PFIX) {
    PF <- DX %*% MO$output$BUGSoutput$mean$B    # linear predictor
    PF <- cbind(rep(1,S), exp(PF))              # OdSTUD 
    PF <- cbind(STUD[, ID],                     # added ID variables 6 Nov 2019
                as.data.frame(PF/rowSums(PF)))        # Proportions
    # names(PF) <- c(ID, paste(TC, LA, RE, 'F', sep='.'))             # Put names on variables / added ID names & changed output names 6 Nov 2019
    names(PF) <- c(ID, TC)  # changed by DP (06/03/20)
    
  } else PF <- NULL
  
  # PREDICT True causes distribution WITH RANDOM EFFECTS
  if (PRAN) {
    PR <- MO$output$BUGSoutput$mean$re[STUD$rG, ]   # Prepare RE matrix
    MR <- colMeans(MO$output$BUGSoutput$mean$re, na.rm = T)    # unweighted means of RE
    if (REME == 0) MR <- rep(0, length(MR))
    PR[is.na(PR)] <- MR[ceiling(which(is.na(PR))/nrow(PR))] # Put MR if no RE for study
    PR <- DX %*% MO$output$BUGSoutput$mean$B + PR # linear PRedictor
    PR <- cbind(rep(1,S), exp(PR))   # OdSTUD 
    #PR <- cbind(STUD[, paste(groupcat)],
    PR <- cbind(STUD[, c(paste(groupcat), ID)],      # added ID variables 6 Nov 2019
                as.data.frame(PR/rowSums(PR))) # PRoportions
    # names(PR) <- c(groupcat, ID, paste(TC, LA, RE, 'R', sep='.'))  # Put names on variables / added ID variables 6 Nov 2019
    names(PR) <- c(groupcat, ID, TC)  # changed by DP (06/03/20)
  } else PR <- NULL
  
  ## PREDICT OBSERVED DEATHS (If data available)
  if(!is.null(DEAT)){
    # Count number of appearences of each sid
    nsid <- as.vector(table(DEAT$sid))
    # Vector with number of cases in each sid
    nn <- as.vector(tapply(DEAT$n, DEAT$sid, sum))
    nn <- rep(nn, nsid)
    if (PFIX) {
      # Build matrices with replacited PF values
      PFmat <- PF[rep(seq_len(nrow(PF)), times = nsid), ]
      PFmat <- PFmat[order(PFmat$sid), ]
      # Fixed random effects
      DEAT[, paste0("pf",PSUF)] <- nn * apply(DEAT[, TC] * PFmat[, TC], 1, sum)
    } 
    if (PRAN) {
      # Build matrices with replacited PR values
      PRmat <- PR[rep(seq_len(nrow(PR)), times = nsid), ]
      PRmat <- PRmat[order(PRmat$sid), ]
      # Random effects predictions
      DEAT[, paste0("pr",PSUF)] <- nn * apply(DEAT[, TC] * PRmat[, TC], 1, sum)
    }
    
    # Was study in estimation sample?
    DEAT[, paste0("e.samp",PSUF)] <- DEAT$sid %in% MO$Studies$sid
    DEAT$p.name <- MO$param$name
    for(v in ADPA) DEAT[,paste0("p.",v)] <- MO$param[[v]]
  }
  
  return(list(Para_pre = match.call(), Para_est = MO$param,
              Pfix = PF, Pran = PR, Pdat = DEAT))
}


###
###  END OF FUNCTION
###
################################################ #





################################################ #
####
####   PREDICTION WITH CI by PANCHO
####  


f.ci <- function(out, dat.pred, id=c("isocode","year"), interval=0.95, PFIX=T, PRAN=F, nset=NULL, perProp=c(0, 0)){
  # out       BUGS object with all the sets of beta estimates from the MCMC
  # dat.pred  data set with covariates to be used in the prediction
  # id        variables in dat.pred that uniquely identify observations
  # interval  Confidence interval. Default 95%
  # PFIX      Fixed effects TRUE/FALSE
  # PRAN      Random effects TRUE/FALSE
  # nset      Number of sets to estimate confidence intervals (optional)
  # perProp   Period Proportions ('early' and 'late')
  
  # Check the vector of period proportions is correct
  if (any(perProp < 0) | length(perProp) != 2) {
    stop('The vector of period proportions is not correct')
  }
  # Normalize Period Proprtions
  if (!sum(perProp) %in% c(0, 1)) {
    perPropNew <- round(perProp / sum(perProp), 2)
    warning(paste0("The vector of period proportions (", perProp[1], ", ", perProp[2], ") has been normalized to (", perPropNew[1], ", ", perPropNew[2], ")"))
    perProp <- perProp / sum(perProp)
  }
  
  # Recover simulation parameters from BUGS output
  niter <- out$output$n.iter
  nchain <- out$output$BUGSoutput$n.chains
  burnin <- out$output$BUGSoutput$n.burnin
  thinin <- out$output$BUGSoutput$n.thin
  ncovs <- nrow(out$output$BUGSoutput$mean$B)
  ncause <- ncol(out$output$BUGSoutput$mean$B)
  nre <- nrow(out$output$BUGSoutput$mean$re)
  
  # Parameters to keep from the simulations
  keep <- rev(seq(niter, burnin+1, -thinin))  # DP: no clear why use this thinin here!
  # Limit the number of sets considered, if 'nset' informed
  if (!is.null(nset)) {
    if (nset >= nchain & length(keep) > nset/nchain) {
      keep <- round(rev(seq(niter, burnin+1, 
                            length.out = round(nset/nchain))))
    }
  }
  
  # Select coefficients from BUGS output (array)
  B <- which(substr(dimnames(out$output$BUGSoutput$sims.array)[[3]], 1, 2) == 'B[')
  B <- out$output$BUGSoutput$sims.array[,,B]
  
  # Select random effects from BUGS output (array)
  if (PRAN) {
    re <- which(substr(dimnames(out$output$BUGSoutput$sims.array)[[3]], 1, 3) == 're[')
    re <- out$output$BUGSoutput$sims.array[,,re]
  }
  
  
  # Auxilary list with some objects of the model output
  out1 <- list()
  out1$param <- out$param
  out1$Vars <- out$Vars
  out1$Studies <- out$Studies
  
  
  # Array to store estimated proportions
  estProp <- array(NA, dim = c(nrow(dat.pred), 
                               length(out$param$VDT), 
                               length(keep)*nchain))
  # Variable to cotrol the 3rd dimension of the array.
  k <- 1    
  
  # PREDICTION
  # for (i in keep) { #########AM: changed 29/7/19 because keep does not index properly when there is burnin
  for (i in 1:length(keep)) {
    for (j in 1:nchain) {
      
      # Trick the function: use the selected set of betas as 'fake mean'
      out1$output$BUGSoutput$mean$B <- matrix(B[i, j, ], nrow = ncovs, ncol = ncause)
      # Trick the function: use the selected set of RE as 'fake mean'
      if (PRAN) {
        out1$output$BUGSoutput$mean$re <- matrix(re[i, j, ], nrow = nre, ncol = ncause)
      } else out1$output$BUGSoutput$mean$re <- NULL
      
      # Call f.p2() function
      if (sum(perProp) == 0) {
        outPred <- f.p2(MO = out1, STUD = dat.pred, ID = id, PFIX = PFIX, PRAN = PRAN)
      } else {
        # Early
        outPred1 <- f.p2(MO = out1, STUD = dat.pred, ID = id, PFIX = PFIX, 
                         PRAN = PRAN, PERIOD = 'early')
        # Late and total
        outPred <- f.p2(MO = out1, STUD = dat.pred, ID = id, PFIX = PFIX, 
                        PRAN = PRAN, PERIOD = 'late')
        if (!is.null(outPred$Pfix)) {
          outPred$Pfix <- outPred1$Pfix * perProp[1] + outPred$Pfix * perProp[2]
        }
        if (!is.null(outPred$Pran)) {
          outPred$Pran <- outPred1$Pran * perProp[1] + outPred$Pran * perProp[2]
        }
        if (!is.null(outPred$Pdat)) {
          outPred$Pdat <- outPred1$Pdat * perProp[1] + outPred$Pdat * perProp[2]
        }
        rm(outPred1)
      }
      
      # Save predictions from current set of estimates
      if (!is.null(outPred$Pran)) {
        #estProp[,,k] <- as.matrix(outPred$Pran)
        estProp[,,k] <- as.matrix(outPred$Pran[,grep(paste(out$param$VDT,collapse="|"), 
                                                     names(outPred$Pran))]) ##27/8/2019 remove random effect variable from estimate matrix. Note this also overrides fixed effects estimates if both are requested.
        ##06/11/2019 remove non-estimate variables
      } else estProp[,,k] <- as.matrix(outPred$Pfix[,grep(paste(out$param$VDT,collapse="|"), 
                                                          names(outPred$Pfix))]) ##06/11/2019 remove non-estimate variables
      k <- k + 1
      
      # Capture identification indices (note if future functionality allows both PFIX and PRAN to be
      # run in one call, this will have to be updated) 6 Nov 2019
      if (i==1 & j==1 & PRAN==TRUE) masterID  <- outPred$Pran[,id]
      if (i==1 & j==1 & PRAN==TRUE) cod.names <- names(outPred$Pran)[grep(paste(out$param$VDT,collapse="|"), names(outPred$Pran))]
      if (i==1 & j==1 & PFIX==TRUE) masterID  <- outPred$Pfix[,id]
      if (i==1 & j==1 & PFIX==TRUE) cod.names <- names(outPred$Pfix)[grep(paste(out$param$VDT,collapse="|"), names(outPred$Pfix))]
      
      rm(outPred)
    }
  }
  
  # Upper and lower bounds of the confidence interval
  if (interval >= 1 | interval <= 0) interval <- 0.95
  CI <- c((1 - interval)/2, (1 + interval)/2)
  
  # Mean, meadian and quantiles (addeed MasterID & changed names 6 Nov 2019)
  Q1 <- as.data.frame(cbind(masterID, apply(estProp, c(1, 2), quantile, CI[1])))
  colnames(Q1) <- c(names(masterID), cod.names)
  Qmean <- as.data.frame(cbind(masterID, apply(estProp, c(1, 2), mean)))
  colnames(Qmean) <- c(names(masterID), cod.names)
  Q2 <- as.data.frame(cbind(masterID, apply(estProp, c(1, 2), median)))
  colnames(Q2) <- c(names(masterID), cod.names)
  Q3 <- as.data.frame(cbind(masterID, apply(estProp, c(1, 2), quantile, CI[2])))
  colnames(Q3) <- c(names(masterID), cod.names)
  
  # Output
  return(list(Q1 = Q1, mean = Qmean, median = Q2, Q3 = Q3, 
              estProp = estProp))
  
}


###
###  END OF FUNCTION
###
################################################ #




################################################ #
####
####   Repeated predictions with MCMC object by DAVID
####  


f.pr <- function(MO, STUP, NP=500, IDV=c("isocode","year"), PPR=c(0, 0)){
  # MO    BUGS object with all the sets of beta estimates from the MCMC
  # STUP  data set with covariates to be used in the prediction
  # IDV   variables in STUP that uniquely identify observations
  # NP    Number of sets to estimate confidence intervals
  # PPR   Period Proportions ('early' and 'late')
  
  S   <- dim(STUP)[1]    # number of studies to predict  
  # extract relevant components from the prediction object
  MA  <- MO$Vars # array of means of covariates
  VDT <- MO$param$VDT
  SA  <- MO$output$BUGSoutput$sims.array  # array of estimations
  # Recover simulation parameters from BUGS output
  VXF <- MA$xvar
  K   <- length(VXF) # number of explanatory variables (excluding constant)
  C   <- length(VDT) # number of undelying causes of death
  I   <- dim(SA)[1]  # number of iterations
  H   <- dim(SA)[2] # number of chains
  cB  <- which(substr(dimnames(SA)[[3]], 1, 2) == 'B[') # beta coefficients
  cR  <- which(substr(dimnames(SA)[[3]], 1, 3) == 're[') # random effects
  R   <- length(cR)/(C-1) # numer of random effects
  print(paste("K=",K,", C=",C, ", I=",I, ", H=",H, ", R=",R))
  
  # Check the vector of period proportions is correct
  if (any(PPR < 0) | length(PPR) != 2) {
    stop('The vector of period proportions is not correct')
  }
  # Normalize Period Proportions
  if (!sum(PPR) %in% c(0, 1)) {
    PPR <- round(PPR / sum(PPR), 2)
    warning(paste0("The vector of period proportions has been normalized to (", PPR[1], ", ", PPR[2],")"))
  }
  
  ## Prepare predicion data
  #  generate unique identifier for each line
  STUP$UID <- STUP[,IDV[1]]
  if(length(IDV)>1) for(i in IDV[-1]) STUP$UID <- paste(STUP$UID, STUP[,i],sep=".")
  # Prepare early and late variables
  STUP$per.early <- 0
  STUP$per.late  <- 0
  # prepare raw variables dataset from prediction sample
  DX    <- cbind(rep(1,S), as.matrix(STUP[,VXF]))
  meMat <- matrix(MA$mean, nrow=S, ncol=K, byrow=T)
  sdMat <- matrix(MA$sd, nrow=S, ncol=K, byrow=T)
  DX[,-1] <- (DX[,-1] - meMat)/sdMat
  rm(meMat, sdMat)
  # create datasets for early and late
  DXl <- DX
  DXl$per.late <- 1
  DXe <- DX
  DXe$per.early <- 1
  
  # Samples to keep from the simulations
  rc <- sample(H, NP, replace=T)  # random selection of chains
  ri <- sample(c(round(I/2):I), NP, replace=F) # random iterations
  
  # Array to store estimated proportions
  EM <- array(NA, dim = c(S, C, R+1, NP))
  dimnames(EM) <- list(STUP$UID, VDT, c(0:R), paste(rc,ri,sep="."))
  
  # PREDICTION
  for (i in 1:NP){
    # Matrix of coefficients
    BM <- cbind(rep(0,K+1), matrix(SA[ri[i],rc[i],cB], nrow=(K+1) )) #betas
    RM <- cbind(rep(0,R), matrix(SA[ri[i],rc[i],cR], nrow=R )) #rnef
    # Calculate logodds with fixed effects
    if(sum(PPR)==0){
      TM <- EM[,,,i]
      TM[,,1] <- DX %*% BM # logodds with fixed effects
      for(j in 1:R) TM[,,1+j] <- TM[,,1] + matrix(rep(RM[j,],each=S), nrow=S)
      EM[,,,i] <- aperm(apply(TM, c(1,3), function(x) exp(x)/sum(exp(x))), perm=c(2,1,3))
      apply(TM, c(1,3), sum)
    }else{
      TMe <- EM[,,,i]
      TMl <- EM[,,,i]
      TMe[,,1] <- DXe %*% BM # logodds with fixed effects
      TMl[,,1] <- DXl %*% BM # logodds with fixed effects
      for(j in 1:R){
        TMe[,,1+j] <- TMe[,,1] + matrix(rep(RM[j,],each=S), nrow=S)
        TMl[,,1+j] <- TMl[,,1] + matrix(rep(RM[j,],each=S), nrow=S)
      }
      TMe <- aperm(apply(TMe, c(1,3), function(x) exp(x)/sum(exp(x))), perm=c(2,1,3))
      TMl <- aperm(apply(TMl, c(1,3), function(x) exp(x)/sum(exp(x))), perm=c(2,1,3))
      EM[,,,i] <- TMe*PPR[1] + TMl*PPR[2]
    }   
  }
  return(list(VX=MO$Vars, RE=MO$Studies, Summary=MO$output$BUGSoutput$summary, EM=EM))
}


###
###  END OF FUNCTION
###
################################################ #





############################################### #
####
####   CALCULATE CIs from PREDICTIONS by David
####  

f.pci <- function(iso, PM, RT, LA, CI=95){
  # iso  country ISO code
  # PM   predicted object from function f.pc2()
  # RT   average across these random effects predictions (0 means fixed)
  # LA   label for the combination fo random effects
  # CI   Credibility interval
  q <- c((100-CI)/200, (100+CI)/200)
  CY <- grep(iso, dimnames(PM$EM)[[1]], value=T)
  for(cy in CY){
    DD <- PM$EM[cy,,RT+1,]
    pp <- as.data.frame(t(rbind(me=apply(DD, 1, mean), sd=apply(DD, 1, sd), qu=apply(DD, 1, quantile,q))))
    pp$cau  <- rownames(pp)
    pp$iso  <- iso
    pp$year <- as.numeric(gsub(paste0(iso,"."),"", cy))
    if(cy == CY[1]) PP <- pp else PP <- rbind(PP,pp)
  }
  PP$lab <- LA
  return(PP)
}

###
###  END OF FUNCTION
###
################################################ #




################################################ #
####
####   PLOT PREDICTIONS WITH CI by David
####  

f.pp <- function(iso, PM, LA=c("fixed","national","local","global"), ys=c(2000:2015), CA=c("intrapartum", "preterm", "congenital", "sepsis", "pneumonia", "diarrhoea", "tetanus", "other"), CO=c("black", "red", "green", "blue", "aquamarine2", "brown", "darkorchid", "grey"), SH=c(15,17,18,19), YM=NULL, TI=NULL){
  # iso country ISO code
  # PM  matrix of oredictions
  # ys  set of years for prediction
  # LA  label ("fixed", "National", "Local", "Global")
  # CA  Only for these causes of death
  # CO  colours for the causes of death
  # SH  Shape of dots for different studies
  # YM  Maximum Y-value in the plot
  # TI  name of country for title
  if(is.null(TI)) TI <- iso
  PP <- PM[(PM$iso==iso) & (PM$year %in% ys) & (PM$lab %in% LA), ]
  PP$col <- CO[as.numeric(PP$cau)]
  PP$pch <- SH[match(PP$lab, LA)]
  nc     <- length(levels(PP$cau))
  PP     <- PP[PP$cau %in% CA,]
  M <- table(PP$year)[1] # number of points per year
  PP <- PP[order(PP$year, PP$cau, PP$lab),]
  PP$csy <- unlist(by(rep(1,nrow(PP)), PP$year, cumsum))
  PP$xl  <- (PP$year - 0.5) + (PP$csy/(M+1))
  if(is.null(YM)) YM <- max(PP[,4])*1.05
  #par(mar=c(4,4,6,1))
  plot(PP$xl, PP$me, ylim=c(0,YM), pch=PP$pch, las=1, col=PP$col, xlab="Year", ylab="Proportion of deaths", main=paste0(TI, " - Mean [",names(PM)[3]," , ",names(PM)[4],"]"))
  arrows(PP$xl, PP[,3], PP$xl, PP[,4], length=0.05, angle=90, code=3, col=PP$col)
  CC <- unique(PP[,c("cau","lab","col","pch")])
  CC$cl <- paste(CC$cau,substr(CC$lab,1,3),sep=".")
  legend(2000, YM, pch=CC$pch, col=CC$col, legend=CC$cl, ncol=length(LA), yjust=0, xpd=T, cex = 0.7)
  return(PP)
}

###
###  END OF FUNCTION
###
############################################### #















