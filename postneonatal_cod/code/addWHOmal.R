
addWHOmal<- function(dat) {
  # dat is one MCMC iteration of predicted proportions for national data
  #i = simulation number

  ##Envelopes
  #  temp<- array2df(results.jtc[i,,])
  #  temp$iso3 <- as.character(temp$d2)
  #  temp$year <- as.numeric(as.character(temp$d1)) - 0.5
  #  temp$u5mr  <- temp[,1]
  #  summary(temp)
  #  u5mr_sample<- temp[which(temp$iso3 %in% data.predict$isocode
  #                               & temp$year %in% data.predict$year
  #                               ), c("iso3","year","u5mr")]
  #  dim(u5mr_sample)
  #  summary(u5mr_sample)

  #triangular distribution
  u_mal<- runif(1,0,1)
  malunc$malaria_sample <- ifelse( u_mal< (malunc$est - malunc$lower) /(malunc$upper - malunc$lower),
                                   malunc$lower +sqrt(u_mal*(malunc$upper - malunc$lower)*(malunc$est - malunc$lower)),
                                   malunc$upper-sqrt((1-u_mal)*(malunc$upper-malunc$lower)*(malunc$upper-malunc$est))  )
  malunc$malaria_sample<- ifelse(malunc$upper == 0, 0, malunc$malaria_sample)
  malunc$malaria_sample<- ifelse(malunc$lower==malunc$upper, malunc$lower, malunc$malaria_sample)

  temp<- dat * dp.nat$envelope_aids_measles_free
  temp2<- merge(temp, malunc[which(malunc$malwho==1),c("malaria_sample","malwho")], by="row.names", all.x=T)
  #print(colnames(dat))
  #print(head(temp))
  #print(head(dp.nat))
  #print(colnames(temp2))
  temp2$envelope_aids_measles_free<- rowSums(temp2[,vdr])

  temp2$malaria_sample<- ifelse(is.na(temp2$malaria_sample), temp2$malaria, temp2$malaria_sample)
  #if malaria sample is greater than 70% of envelope, use MCEE estimate
  temp2$malaria_sample<- ifelse(temp2$malaria_sample > temp2$envelope_aids_measles_free*0.7,
                                temp2$malaria, temp2$malaria_sample)
  for(C in vdr[which(vdr != "malaria")]) {
    temp2[,C]<- ifelse(temp2$envelope_aids_measles_free>0, temp2[,C] * (temp2$envelope_aids_measles_free - temp2$malaria_sample) /
                         (temp2$envelope_aids_measles_free - temp2$malaria), temp2[,C])
  }
  temp2$malaria<- temp2$malaria_sample
  rownames(temp2)<- temp2$Row.names
  #cat("\n\n A\n")
  #print(dim(temp2[,vdr]/rowSums(temp2[,vdr])))
  #          print(  dim(dat))
  out<- temp2[,vdr] /rowSums(temp2[,vdr])
  return(out)
  #p21$estR.mal.N[,,i]<- temp2[,vdr]
}
