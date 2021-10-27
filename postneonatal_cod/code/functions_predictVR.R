
##################################################Z
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
#################################################



##################################################Z
####
####   PREDICTION WITH CI FROM DAVID
####


f.ci2 <- function(out, dat.pred, nset=500, id=c("isocode","year"),
                  interval=0.95, PFIX=T, PRAN=F, perProp=c(0, 0)){
  # out       BUGS object with all the sets of beta estimates from the MCMC
  # dat.pred  data set with covariates to be used in the prediction
  # id        variables in dat.pred that uniquely identify observations
  # interval  Confidence interval. Default 95%
  # PFIX      Fixed effects TRUE/FALSE
  # PRAN      Random effects TRUE/FALSE
  # nset      Number of sets to estimate confidence intervals
  # perProp   Period Proportions ('early' and 'late')

  # generate unique identifier for each line
  dat.pred$UID <- dat.pred[,id[1]]
  if(length(id)>1) for(i in id[-1]) dat.pred$UID <- paste(dat.pred$UID, dat.pred[,i],sep=".")

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
  niter  <- dim(out$output$BUGSoutput$sims.array)[[1]] #out$output$BUGSoutput$n.iter
  nchain <- out$output$BUGSoutput$n.chains
  ncovs  <- nrow(out$output$BUGSoutput$mean$B)
  ncause <- ncol(out$output$BUGSoutput$mean$B)
  nre    <- nrow(out$output$BUGSoutput$mean$re)
  vdt    <- out$param$VDT

  # Samples to keep from the simulations
  cat("\nchain ", nchain,"\n")
  nrc <- sample(nchain, nset, replace=T)  # random selection of chains
  nri <- sample(c(round(niter/2):niter), nset, replace=F) # random selection of iterations

  # Select coefficients from BUGS output (array)
  B <- which(substr(dimnames(out$output$BUGSoutput$sims.array)[[3]], 1, 2) == 'B[')
  # Select random effects from BUGS output (array)
  if(PRAN==T) re <- which(substr(dimnames(out$output$BUGSoutput$sims.array)[[3]], 1, 3) == 're[')

  # Auxilary list with some objects of the model output
  out1 <- list()
  out1$param   <- out$param
  out1$Vars    <- out$Vars
  out1$Studies <- out$Studies

  # Arrays to store estimated proportions
  #cat("\nestF\n")
  estF <- array(NA, dim = c(nrow(dat.pred), length(vdt), nset))
  #print(dim(estF))
  #print(summary(nrc))
  #print(summary(nri))
  dimnames(estF) <- list(dat.pred$UID, vdt, paste(nrc,nri,sep="."))

  if(PRAN==T) estR <- estF else estR <- NA

  # PREDICTION
  # for (i in keep) { #########AM: changed 29/7/19 because keep does not index properly when there is burnin
  for (i in 1:nset){
      # Trick the function: use the selected set of betas as 'fake mean'
      out1$output$BUGSoutput$mean$B <- matrix(out$output$BUGSoutput$sims.array[nri[i],nrc[i],B], nrow = ncovs, ncol = ncause)
      # Trick the function: use the selected set of RE as 'fake mean'
      if (PRAN==T) {
        out1$output$BUGSoutput$mean$re <- matrix(out$output$BUGSoutput$sims.array[nri[i],nrc[i],re], nrow = nre, ncol = ncause)
      } else out1$output$BUGSoutput$mean$re <- NULL
      # Call f.p2() function
      if (sum(perProp) == 0) {
   #     cat("\n!!!\n")
        outPred <- f.p2(MO=out1, STUD = dat.pred, ID=id, PFIX=PFIX, PRAN=PRAN, REME=1)
      } else {
        # Early
        outPred1 <- f.p2(MO=out1, STUD=dat.pred, ID=id, PFIX=PFIX, PRAN=PRAN, PERIOD='early')
        # Late and total
        outPred <- f.p2(MO=out1, STUD=dat.pred, ID=id, PFIX=PFIX, PRAN=PRAN, PERIOD='late')
        if (!is.null(outPred$Pfix)) outPred$Pfix <- outPred1$Pfix*perProp[1] + outPred$Pfix*perProp[2]
        if (!is.null(outPred$Pran)) outPred$Pran <- outPred1$Pran*perProp[1] + outPred$Pran*perProp[2]
        if (!is.null(outPred$Pdat)) outPred$Pdat <- outPred1$Pdat*perProp[1] + outPred$Pdat*perProp[2]
        rm(outPred1)
      }
      # Save predictions from current set of estimates
      if (!is.null(outPred$Pran)) estR[,,i] <- as.matrix(outPred$Pran[,vdt]) ##06/11/2019 remove non-estimate variables
      if (!is.null(outPred$Pfix)) estF[,,i] <- as.matrix(outPred$Pfix[,vdt]) ##06/11/2019 remove non-estimate variables
      rm(outPred)
    }

  # Output
  return(list(estF=estF, estR=estR))
}


###
###  END OF FUNCTION
###
#################################################




