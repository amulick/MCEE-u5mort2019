
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
  ##  -------------------------
  ##  STUD = dataset of studies
  ##  DEAT = dataset of reported causes of death by study with missclassification matrix
  ##  VDT  = names of true causes of death (string vector)
  ##  VXF  = names of fixed effects variables (string vector)
  ##  VXR  = name of random effects variable (string vector)
  ##  --------------
  ##  MODL = txt file with BUGS model
  ##  LAMB = Lambda value for Poisson
  ##  RSDL = limit for the random effects SD
  ##  QUAD = NULL if no quadratic terms desired (default); !is.null(QUAD) adds quadratic terms
  ##  --------------
  ##  PARA = 0 run CHAINS in serial, PARA!=0 run CHAINS in parallel
  ##  CHAS = number of mcmc SAMPLES
  ##  ITER = number of draws used for burn-in (currently arbitrarily set to 500)
  ##  BURN = number of chains
  ##  THIN = thin factor 
  ##  SEED = random seed for jags.parallel (for jags use set.seed() first)
  ##  --------------
  ##  PRIN = Print stage of computation
  ##  SAMC = Save MCMC samples? TRUE/FALSE
  ##  SAVX = T, Save explanatory variables
  ##  SADE = T, Save matrix of deaths
  ##  NAME = Optional name for the model
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
#################################################




##################################################Z
####
####   Prediction from ONE model estimation in other dataset
####  

f.p1 <- function(MO, STUD, DEAT, ADPA=NULL, PFIX=T, PRAN=T, REME=0, PSUF=NULL){
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

  # Add random effect group to the prediction data
  STUD <- merge(STUD, unique(MO$Studies[,c(2,3)]), by=names(MO$Studies)[2], all.x=T, all.y=F)
  # prepare raw variables dataset from prediction sample
  DX <- cbind(rep(1,S), as.matrix(STUD[,MO$Vars$xvar]))
  # normalise prediction variables using means and SSTUD of estimation dataset
  meanMat <- matrix(MO$Vars$mean, nrow = nrow(DX), ncol = K, byrow = T)
  sdMat <- matrix(MO$Vars$sd, nrow = nrow(DX), ncol = K, byrow = T)
  DX[,-1] <- (DX[,-1] - meanMat)/sdMat
  rm(meanMat, sdMat)

  # Prepare average andom effects
  PR <- MO$output$BUGSoutput$mean$re[STUD$rG, ]   # Prepare RE matrix
  MR <- colMeans(MO$output$BUGSoutput$mean$re, na.rm = T)    # unweighted means of RE
  #print("HEY!")
  #print(dim(PR))
  #print(head(PR))
  
  if(REME == 0) MR <- rep(0, length(MR))
  PR[is.na(PR)] <- MR[ceiling(which(is.na(PR))/nrow(PR))] # Put MR if no RE for study
  #print("HEY 2!")
  #print(dim(PR))
  #print(head(PR))
  
  # PREDICTION True causes distribution WITH FIXED EFFECTS
  if(PFIX) {
    PF <- DX %*% MO$output$BUGSoutput$mean$B    # linear predictor
    PF <- cbind(rep(1,S), exp(PF))  # OdSTUD
    PF <- cbind(sid=STUD$sid, as.data.frame(PF/rowSums(PF))) # Proportions
    names(PF)[-1] <- TC             # Put names on variables
  }
  #print(dim(DX))
  #print(dim(PF))

  # PREDICT True causes distribution WITH RANDOM EFFECTS
  if(PRAN) {
    PR <- DX %*% MO$output$BUGSoutput$mean$B + PR # linear PRedictor
    PR <- cbind(rep(1,S), exp(PR))   # OdSTUD
    PR <- cbind(sid=STUD$sid, as.data.frame(PR/rowSums(PR))) # PRoportions
    names(PR)[-1] <- TC  # Put names on variables
  }
  #print("HEY 3!")
  ## PREDICT OBSERVED DEATHS (If data available)
  if(!is.null(DEAT)){
    # Number of appearences of each sid
    nsid <- as.vector(table(DEAT$sid))
    # Vector with number of cases in each sid
    nn <- as.vector(tapply(DEAT$n, DEAT$sid, sum))
    nn <- rep(nn, nsid)
    # fixed effects predictions
    if (PFIX) {
      PFmat <- PF[match(unique(DEAT$sid), PF$sid), ]
      PFmat <- PFmat[rep(seq_len(nrow(PFmat)), times = nsid), ]
      DEAT[, paste0("pf",PSUF)] <- nn * rowSums(DEAT[, TC] * PFmat[, TC])
    }
    # random effects predictions
    if (PRAN) {
      # PRmat <- PR[rep(seq_len(nrow(PR)), times = nsid), ]
      # PRmat <- PRmat[order(PRmat$sid), ]
      PRmat <- PR[match(unique(DEAT$sid), PR$sid), ]
      PRmat <- PRmat[rep(seq_len(nrow(PRmat)), times = nsid), ]
      DEAT[, paste0("pr",PSUF)] <- nn * rowSums(DEAT[, TC] * PRmat[, TC])
    }
    # Was study in estimation sample?
    DEAT[, paste0("e.samp",PSUF)] <- DEAT$sid %in% MO$Studies$sid
    # Add names
    DEAT$p.name <- MO$param$name
    for(v in ADPA) DEAT[,paste0("p.",v)] <- MO$param[[v]]
  }

  # Output
  rtn <- list()
  rtn[["Para_pre"]] <- match.call()
  rtn[["Para_est"]] <- MO$param
  if(!is.null(PF)) rtn[["Pfix"]] <- PF
  if(!is.null(PR)) rtn[["Pran"]] <- PR
  if(!is.null(PF)) rtn[["Pdat"]] <- DEAT
  return(rtn)
}

###
###  END OF FUNCTION
###
#################################################Z



#################################################
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
  ##  ----------------------------------------
  ##  STUD = dataset of studies
  ##  DEAT = dataset of reported causes of death by study with missclassification matrix
  ##  VDT  = names of true causes of death (string vector)
  ##  VXF  = names of fixed effects variables (string vector)
  ##  VXR  = name of random effects variable (string vector)
  ##  --------------
  ##  LAMB = Lambda value for Poisson
  ##  RSDL = limit for the random effects SD
  ##  QUAD = NULL if no quadratic terms desired (default); !is.null(QUAD) adds quadratic terms
  ##  --------------
  ##  MODL = txt file with BUGS model
  ##  PARA = 0 run CHAINS in serial, PARA!=0 run CHAINS in parallel
  ##  CHAS = number of mcmc SAMPLES
  ##  ITER = number of draws used for burn-in (currently arbitrarily set to 500)
  ##  BURN = number of chains
  ##  THIN = thin factor 
  ##  ---------------
  ##  PRIN = Print stage of computation
  ##  SAMC = Save JAGS object for MCMC diagnostics? TRUE/FALSE
  ##  SAVX = T, Save explanatory variables
  ##  SADE = T, Save matrix of deaths
  ##  CORE = Number of cores to use for CV samples (>1 = parallelise CV partitions (not chains in mcmc (default)), 1=parallel chains (not CV partitions), 0= do not parallelise anything)
  ##  ---------------
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
#################################################Z




##################################################Z
####
####   Prediction from a set of CV results
####  

# CVS <- vanams; STUD <- studies; DEAT <- deaths
# ADPA <- c("LAMB","RSDL"); REME <- 0; PFIX <- T; PRAN <- T

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
    cat("cv ::: ",cv,"/n")  
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
  tt <- aggregate(PT[,c("All")], lvt, length) # count number of studies in strata
  names(tt) <- c(vt, "N")
  if(is.null(PT$pf)==F){
    PT$pf1 <- (PT$pf - PT$n)^2
    PT$pf2 <- (PT$pf1)/PT$pf
    trf <- aggregate(PT[,c("pf2")], lvt, function(x) sqrt(sum(x))) #root mean square error
    names(trf) <- c(vt,"pf2")
    for(c in c("pf2")) trf[,c] <- formatC(trf[,c], 1, format="f")
    tt <- merge(tt,trf, by=vt)
  }
  
  if(is.null(PT$pr)==F){
    PT$pr1 <- (PT$pr - PT$n)^2
    PT$pr2 <- (PT$pr1)/PT$pr
    trr <- aggregate(PT[,c("pr2")], lvt, function(x) sqrt(sum(x))) #root mean square error
    names(trr) <- c(vt,"pr2")
    for(c in c("pr2")) trr[,c] <- formatC(trr[,c], 1, format="f")
    tt <- merge(tt,trr, by=vt)
  }
  return(tt)
}

###
###  END OF FUNCTION
###
#################################################Z




