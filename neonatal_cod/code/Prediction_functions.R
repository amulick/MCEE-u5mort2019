

## Functions needed for Prediction, By David Prieto




################################################ #
####
####   Object of arrays of coefficients from MCMC

f.par <- function(MO, NP=500){
  # MO    BUGS object with all the sets of beta estimates from the MCMC
  # NP    Number of sets to estimate confidence intervals
  # This  makes cubes of fixed and random effects coefficients
  # extract relevant components from the prediction object
  VXF  <- MO$Vars$xvar # array of means of covariates
  VDT <- MO$param$VDT
  SA  <- MO$output$BUGSoutput$sims.array  # array of estimations
  # Recover simulation parameters from BUGS output
  K   <- length(VXF) # number of explanatory variables (excluding constant)
  C   <- length(VDT) # number of underlying causes of death
  I   <- dim(SA)[1]  # number of iterations
  H   <- dim(SA)[2] # number of chains
  cB  <- which(substr(dimnames(SA)[[3]], 1, 2) == 'B[') # beta coefficients
  cR  <- which(substr(dimnames(SA)[[3]], 1, 3) == 're[') # random effects
  R   <- length(cR)/(C-1) # number of random effects
  print(paste("K=",K,", C=",C, ", I=",I, ", H=",H, ", R=",R))
  # Samples to keep from the simulations
  rc <- sample(H, NP, replace=T)  # random selection of chains
  ri <- sample(I, NP, replace=F) # random iterations
  # Prepare arrays
  BM <- array(NA, dim = c(K+1, C, NP))  # array of fixed betas
  dimnames(BM) <- list(c("constant",VXF), VDT, paste(rc,ri,sep="."))
  RM <- array(NA, dim = c(R, C, NP))  # array of random effects
  dimnames(RM) <- list(MO$Studies[match(c(1:R), MO$Studies$rG), 2], VDT, paste(rc,ri,sep="."))
  # populate arrays
  for (i in 1:NP){
    # Matrix of coefficients
    BM[,,i] <- cbind(rep(0,K+1), matrix(SA[ri[i],rc[i],cB], nrow=(K+1) )) #betas
    RM[,,i] <- cbind(rep(0,R),   matrix(SA[ri[i],rc[i],cR], nrow=R)) #rnef
  }
  return(list(Model=MO$param$name, VX=MO$Vars, VDT=VDT, RE=MO$Studies, Summary=MO$output$BUGSoutput$summary, BM=BM, RM=RM))
}

###
###  END OF FUNCTION
###
################################################ #


################################################ #
####
####   Predictions with an array of coefficients from function f.par()
####   For early/late just create variables


f.pr2 <- function(PA, PD, RT, ID="id", PE=""){
  # PA   object produced with function f.par()
  # PD   data set with covariates to be used in the prediction
  # ID   variable in PD that uniquely identify observations
  # PE   Period label ("early","late","any")
  # RT  list of applicable random terms for each ID in PD
  # This function makes predictions with fixed effects only and then adds ONE random effect term selected at random among some canditates in RT in each mcmc iteration.
  # generate unique identifier for each line
  ## Prepare prediction data
  VXF <- PA$VX$xvar  # list of prediction covariates
  S   <- dim(PD)[1]
  K   <- dim(PA$BM)[1]-1 # number of covariates excluding constant
  C   <- dim(PA$BM)[2]   # number of causes of death
  N   <- dim(PA$BM)[3]   # number of simulations
  if(sum(!(VXF %in% names(PD)))>0){
    stop(paste('Variables', paste(VXF[!(VXF %in% names(PD))], collapse=" , "), ' not in the data'))
  }
  # prepare raw variables dataset from prediction sample
  for(k in 1:K) PD[,VXF[k]] <- (PD[,VXF[k]] - PA$VX$mean[k])/PA$VX$sd[k]
  DX <- cbind(cons=rep(1,S), as.matrix(PD[,VXF]))
  # Prepare matrix of selected random terms for each ID
  RE <- array(NA, dim = c(S, N))
  #  for(i in 1:S) RE[i,] <- ifelse(length(RT[[i]])>0, sample(RT[[i]], N, replace=T), sample(dimnames(PA$RM)[[1]], N, replace=T))

  #Corrections 8/13/2020
  for(i in 1:S) {
    if(length(RT[[i]]>0)) {
    #Corrections 8/27/2020 - because when RT[[i]] has length 1 and is numeric, the function samples from 1:RT[[i]] rather than just RT[[i]]
      # RE[i,] <- sample(RT[[i]], N, replace=T)
      RE[i,] <- sample(as.character(RT[[i]]), N, replace=T)
    } else {
      RE[i,]  <- sample(dimnames(PA$RM)[[1]], N, replace=T)
    }
  }

  # Prepare matrix to store predictions
  LF <- array(NA, dim = c(S, C, N))  # matrix of fixed Predictions
  dimnames(LF) <- list(PD[,ID], dimnames(PA$BM)[[2]], dimnames(PA$BM)[[3]])
  LR <- LF  # matrix of random Predictions
  # PREDICTION for each iteration
  for (i in 1:N){
    LF[,,i] <- DX %*% PA$BM[,,i] # logodds with fixed effects
    LR[,,i] <- LF[,,i] + PA$RM[RE[,i],,i] # logodds of random effects
  }
  PF <- aperm(apply(LF, c(1,3), function(x) exp(x)/sum(exp(x))), perm=c(2,1,3))
  PR <- aperm(apply(LR, c(1,3), function(x) exp(x)/sum(exp(x))), perm=c(2,1,3))
  return(list(Model=PA$Model, Period=PE, VX=PA$VX, PF=PF, PR=PR) )
}

###
###  END OF FUNCTION
###
################################################ #






############################################## #
###
###   CALCULATE CIs from PREDICTIONS by David
###

 f.pci2 <- function(PM, CI=95){
   # PM   predicted object from function f.pr2() by David
   # CI   Credibility interval
   q <- c((100-CI)/200, (100+CI)/200)
   for(c in 1:(dim(PM$PF)[2])){
     xf <- data.frame(t(apply(PM$PF[,c,], 1, function(x) c(mean(x, na.rm=T), sd(x, na.rm=T), quantile(x, q[1], na.rm=T), quantile(x, q[2], na.rm=T)))))
     xf$cause <- dimnames(PM$PF)[[2]][c]
     xf$type <- "fixed"
     xf$iso  <- sapply(rownames(xf), function(x)  strsplit(x, ".",fixed=T)[[1]][1])
     xf$year <- as.numeric(sapply(rownames(xf), function(x)  strsplit(x, ".",fixed=T)[[1]][2]))
     xr <- data.frame(t(apply(PM$PR[,c,], 1, function(x) c(mean(x, na.rm=T), sd(x, na.rm=T), quantile(x, q[1], na.rm=T), quantile(x, q[2], na.rm=T)))))
     xr$cause <- dimnames(PM$PR)[[2]][c]
     xr$type <- "random"
     xr$iso  <- sapply(rownames(xr), function(x)  strsplit(x, ".",fixed=T)[[1]][1])
     xr$year <- as.numeric(sapply(rownames(xr), function(x)  strsplit(x, ".",fixed=T)[[1]][2]))
     if(c==1) PP <- rbind(xf,xr) else PP <- rbind(PP,xf,xr)
   }
   names(PP)[1:4] <- c("me","sd","lo","hi")
   PP$model <- PM$Model
   PP$period <- PM$Period
   return(PP[,c(7:10,6,5,1:4)])
 }

###
###  END OF FUNCTION
###
################################################ #









################################################ #
####
####   PLOT PREDICTIONS WITH CI by David
####

f.pp2 <- function(ISO, PM, PD="", MO="", LA=c("fixed","random"), YS=c(2000:2019), CA="", YM=NULL, TI=NULL){
  # ISO country iso code
  # PM  matrix of predictions from function f.pci2()
  # YS  set of years for prediction
  # LA  label ("fixed", "National", "Local", "Global")
  # CA  Only for ONE causes of death
  # YM  Maximum Y-value in the plot
  # TI  name of country for title
  print(MO)
  if(mean(MO %in% unique(PM$model))<1) stop(paste('Need to specify ALL the models (MO) you want to plot:', paste(unique(PM$model), collapse=", ")))
  if(!(PD %in% unique(PM$period))) stop(paste('Need to specify ONE Period of Deaths (PD) existing in data:', paste(unique(PM$period), collapse=", ")))
  if(!(CA %in% unique(PM$cause))) stop(paste('Need to specify ONE Cause of Deaths (CA) existing in data:', paste(unique(PM$cause), collapse=", ")))
  if(is.null(TI)) TI <- ISO
  D1 <- filter(PM, iso==ISO, period==PD, cause==CA, model %in% MO, year %in% YS, type %in% LA) %>%  mutate(xp=year+0.2*(model=="VR")+0.4*(model=="VAVR"))
  if(is.null(YM)) YM <- max(D1$hi, na.rm=T)*1.05
  ggplot(D1, aes(x=year, col=model, shape=type)) + ylim(0,YM) +
        geom_point(aes(x=xp, y=me)) +
        geom_errorbar(aes(x=xp, ymin=lo, ymax=hi)) +
        labs(title=paste("Country: ",TI,"\nPeriod",PD," - Cause:",CA), subtitle=paste("Comparison of predictions with models:", paste(MO, collapse=" and ")), y="Estimated proportion of deaths") + theme(legend.position="bottom")
}

###
###  END OF FUNCTION
###
############################################### #




################################################ #
####
####   PLOT comparison of models (6 causes, all years in a country)
####

f.pcm<- function(ISO, PM, PD="", MO="", LA=c("fixed","random"), YS=c(2000:2015), YM=NULL, TI=NULL, CA=c("intrapartum", "preterm", "congenital", "sepsis", "pneumonia", "other")){
  # ISO country iso code
  # PM  matrix of predictions from function f.pci2()
  # PD  Period of death
  # MO  Model (VA, VR)
  # YS  set of years for prediction
  # LA  label ("fixed", "National", "Local", "Global")
  # YM  Maximum Y-value in the plot
  # TI  name of country for title
  # CA list of causes to plot
  if(mean(MO %in% unique(PM$model))<1) stop(paste('Need to specify ALL the models (MO) you want to plot:', paste(unique(PM$model), collapse=", ")))
  if(!(PD %in% unique(PM$period))) stop(paste('Need to specify ONE Period of Deaths (PD) existing in data:', paste(unique(PM$period), collapse=", ")))
  if(is.null(TI)) TI <- ISO
  PM <- filter(PM, iso==ISO, period==PD, model %in% MO, year %in% YS) %>% mutate(xp=year+0.2*(model=="VR")+0.4*(type=="random"))
  lpo <- list()
  for(c in 1:length(CA)){
    D1 <- filter(PM, cause==CA[c])
    if(is.null(YM)) ym <- max(D1$hi, na.rm=T)*1.05 else ym <- YM
    lpo[[c]] <- ggplot(D1, aes(x=year, col=model, shape=type)) + ylim(0,ym) +
      geom_point(aes(x=xp, y=me)) +
      geom_errorbar(aes(x=xp, ymin=lo, ymax=hi, width=0.3), linetype=5) +
      ylab("Estimated proportion of deaths")
  }
  pg <- annotate_figure(ggarrange(plotlist=lpo, nrow=3, ncol=2, labels=CA, common.legend=T, legend="bottom"), top=paste("Country: ",TI,"\nPeriod of deaths:",PD," - Comparison of predictions with models:", paste(MO, collapse=" and ")))
  return(pg)
}

###
###  END OF FUNCTION
###
############################################### #



################################################ #
####
####   PLOT PREDICTIONS WITH CI by David
####

f.pmp <- function(PM, PD="", MO=c("VA","VR"), LA="random", YS=2019, CA="", YM=NULL){
  # PM  matrix of predictions from function f.pci2()
  # PD  Period of death
  # MO  Model (VA, VR)
  # YS  set of years for prediction
  # LA  label ("fixed", "National", "Local", "Global")
  # CA  Only for these causes of death
  # YM  Maximum Y-value in the plot
  # TI  name of country for title
  if(mean(MO %in% unique(PM$model))<1) stop(paste('Need to specify ALL the models (MO) you want to plot:', paste(unique(PM$model), collapse=", ")))
  if(!(PD %in% unique(PM$period))) stop(paste('Need to specify ONE Period of Deaths (PD) existing in data:', paste(unique(PM$period), collapse=", ")))
  if(!(CA %in% unique(PM$cause))) stop(paste('Need to specify ONE Cause of Deaths (CA) existing in data:', paste(unique(PM$cause), collapse=", ")))
  PM <- filter(PM, type==LA, period==PD, cause==CA, model %in% MO, year==YS) %>%  mutate(xp=as.numeric(as.factor(iso))+0.2*(model=="VR"))
  if(is.null(YM)) YM <- max(PM$hi)*1.05
  ggplot(PM, aes(x=me, y=xp, col=model)) + xlim(0,YM) + geom_point() +
    geom_errorbarh(aes(xmin=lo, xmax=hi), linetype=1) +
    scale_y_discrete(name="",limits=levels(as.factor(PM$iso)), labels=levels(as.factor(PM$iso))) +
    labs(title=paste("Period",PD," - Cause:",CA, "- Year:",YS, "- Effects:",LA), x="Estimated proportion of deaths", color="Comparison of predictions with models:") + theme(legend.position="top")
}

###
###  END OF FUNCTION
###
############################################### #












########################################################## #
###
###   OLD FUNCTIONS
###


################################################ #
####
####   Repeated predictions with MCMC object by DAVID
####   For early/late just create variables
#
# f.pr2 <- function(MO, STUP, NP=500, IDV=c("isocode","year"), PD=""){
#   # MO    BUGS object with all the sets of beta estimates from the MCMC
#   # STUP  data set with covariates to be used in the prediction
#   # NP    Number of sets to estimate confidence intervals
#   # IDV   variables in STUP that uniquely identify observations
#   # PD    Period label ("early","late","any")
#   # This function makes predictions with fixed effects only and then adds ONE random effects terms (selected at random). So within each mcmc iteration selected, you have 2 sets of predictions for each country.
#   if(!(PD %in% c("early","late","any"))){
#     stop('Need to specify Period of Deaths parameter as one of PD=c("early","late","any")')
#   }
#   S   <- dim(STUP)[1]    # number of studies to predict
#   # extract relevant components from the prediction object
#   MA  <- MO$Vars # array of means of covariates
#   VDT <- MO$param$VDT
#   SA  <- MO$output$BUGSoutput$sims.array  # array of estimations
#   # Recover simulation parameters from BUGS output
#   VXF <- MA$xvar
#   K   <- length(VXF) # number of explanatory variables (excluding constant)
#   C   <- length(VDT) # number of undelying causes of death
#   I   <- dim(SA)[1]  # number of iterations
#   H   <- dim(SA)[2] # number of chains
#   cB  <- which(substr(dimnames(SA)[[3]], 1, 2) == 'B[') # beta coefficients
#   cR  <- which(substr(dimnames(SA)[[3]], 1, 3) == 're[') # random effects
#   R   <- length(cR)/(C-1) # numer of random effects
#   print(paste("K=",K,", C=",C, ", I=",I, ", H=",H, ", R=",R))
#   ## Prepare prediction data
#   #  generate unique identifier for each line
#   STUP$UID <- STUP[,IDV[1]]
#   if(length(IDV)>1) for(i in IDV[-1]) STUP$UID <- paste(STUP$UID, STUP[,i],sep=".")
#   # Create variables=0 if they are not present in the data but are in the model
#   print(VXF[!(VXF %in% names(STUP))])
#   for(v in VXF) if (!(v %in% names(STUP))) STUP[,v] <- 0
#   print(VXF %in% names(STUP))
#   # prepare raw variables dataset from prediction sample
#   DX    <- cbind(rep(1,S), as.matrix(STUP[,VXF]))
#   meMat <- matrix(MA$mean, nrow=S, ncol=K, byrow=T)
#   sdMat <- matrix(MA$sd, nrow=S, ncol=K, byrow=T)
#   DX[,-1] <- (DX[,-1] - meMat)/sdMat
#   rm(meMat, sdMat)
#   # Samples to keep from the simulations
#   rc <- sample(H, NP, replace=T)  # random selection of chains
#   ri <- sample(c(round(I/4):I), NP, replace=F) # random iterations
#   # Array to store estimated effects
#   FL <- array(NA, dim = c(S, C, NP))  # matrix of fixed log(odds)
#   dimnames(FL) <- list(STUP$UID, VDT, paste(rc,ri,sep="."))
#   RL <- FL  # matrix of random log(odds)
#   # PREDICTION
#   for (i in 1:NP){
#     # Matrix of coefficients
#     BM <- cbind(rep(0,K+1), matrix(SA[ri[i],rc[i],cB], nrow=(K+1) )) #betas
#     RM <- cbind(rep(0,R), matrix(SA[ri[i],rc[i],cR], nrow=R)) #rnef
#     FL[,,i] <- DX %*% BM # logodds with fixed effects
#     RL[,,i] <- FL[,,i] +  RM[sample(c(1:nrow(RM)),S,T),] # logodds with random effects
#   }
#   # Transform in probabilities
#   PF <- aperm(apply(FL, c(1,3), function(x) exp(x)/sum(exp(x))), perm=c(2,1,3))
#   PR <- aperm(apply(RL, c(1,3), function(x) exp(x)/sum(exp(x))), perm=c(2,1,3))
#   return(list(Model=MO$param$name, Period=PD, VX=MO$Vars, RE=MO$Studies, Summary=MO$output$BUGSoutput$summary, PF=PF, PR=PR))
# }
#
###
###  END OF FUNCTION
###
################################################ #



############################################### #
####
####   CALCULATE CIs from PREDICTIONS by David
####
#
# f.pci2 <- function(PM, CI=95){
#   # PM   predicted object from function f.pr2() by David
#   # CI   Credibility interval
#   q <- c((100-CI)/200, (100+CI)/200)
#   for(c in 1:(dim(PM$PF)[2])){
#     xf <- data.frame(t(apply(PM$PF[,c,], 1, function(x) c(mean(x, na.rm=T), sd(x, na.rm=T), quantile(x, q[1], na.rm=T), quantile(x, q[2], na.rm=T)))))
#     xf$cause <- dimnames(PM$PF)[[2]][c]
#     xf$type <- "fixed"
#     xf$iso  <- sapply(rownames(xf), function(x)  strsplit(x, ".",fixed=T)[[1]][1])
#     xf$year <- as.numeric(sapply(rownames(xf), function(x)  strsplit(x, ".",fixed=T)[[1]][2]))
#     xr <- data.frame(t(apply(PM$PR[,c,], 1, function(x) c(mean(x, na.rm=T), sd(x, na.rm=T), quantile(x, q[1], na.rm=T), quantile(x, q[2], na.rm=T)))))
#     xr$cause <- dimnames(PM$PR)[[2]][c]
#     xr$type <- "random"
#     xr$iso  <- sapply(rownames(xr), function(x)  strsplit(x, ".",fixed=T)[[1]][1])
#     xr$year <- as.numeric(sapply(rownames(xr), function(x)  strsplit(x, ".",fixed=T)[[1]][2]))
#     if(c==1) PP <- rbind(xf,xr) else PP <- rbind(PP,xf,xr)
#   }
#   names(PP)[1:4] <- c("me","sd","lo","hi")
#   PP$model <- PM$Model
#   PP$period <- PM$Period
#   return(PP[,c(7:10,6,5,1:4)])
# }

###
###  END OF FUNCTION
###
################################################ #






