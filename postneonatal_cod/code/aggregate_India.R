# Make post hoc adjustment to predicted proportions



aggregate_India<- function(dat) {
# dat is one MCMC iteration of predicted proportions for national data
# columns = causes
# rownames = iso3.year

  causes<- vdr
  #translate %'s to N's
  temp<- dat * data.predict$envelope_aids_measles_free
  dat<- as.data.frame(temp)
  colnames(dat)<- vdr
#  print("PH 1")
#  print(rownames(dat))
  year<- data.predict$year
  iso<- data.predict$isocode
  iso.nat<- substr(iso,1,3)
  #print(table(iso))
  #print(table(iso.nat), useNA="always")
  
  temp<- aggregate(dat, by=list(year, iso.nat), FUN=sum)
  #print(dim(dat))
  #print(dim(temp))
  #print(head(dat))
  rownames(temp)<- paste0(temp$Group.2,".", temp$Group.1)
  #print(head(temp))
  #print(table(temp$Group.1))
  temp2<- temp[,vdr]
  out<- temp2[,vdr] /rowSums(temp2[,vdr])
  return(out)
  
  }
