# Make post hoc adjustment to predicted proportions



posthoc_hib<- function(dat, verbose=FALSE) {
# dat is one MCMC iteration of predicted proportions for national data
# columns = causes
# rownames = iso3.year

  causes<- vdr
  #translate %'s to N's
  temp<- dat * data.predict$envelope_aids_measles_free
  dat<- as.data.frame(temp)
  colnames(dat)<- vdr
#  print("PH 1")
#*******************************************************************************************************
# cause-specific lives saved = cause-specific deaths * serotype% * vaccine coverage * vaccine effectiveness/efficacy
#*******************************************************************************************************

#g pcv_effcov = (sero_pcv7 * cover_pcv7 + sero_pcv10 * cover_pcv10 + sero_pcv13 * cover_pcv13) /100 * 0.8 //0.8 is formulation-independent efficacy for pcv (Lucero 2009)
#label var pcv_effcov "effective coverage for pcv"
  pcv_effcov <- (data.predict$sero_pcv7 * data.predict$cover_pcv7 +
                   data.predict$sero_pcv10 * data.predict$cover_pcv10 +
                   data.predict$sero_pcv13 * data.predict$cover_pcv13) /100 * 0.8
#0.8 is formulation-independent efficacy for pcv (Lucero 2009)

#* take account of the Rotavirus vaccine coverage for diarrhea
delta_dia <- dat$diarrhea * (data.predict$rota/100) * data.predict$dia_rota_effective * data.predict$pdiadeaths_rota
dia_rota  <- dat$diarrhea - delta_dia

#* take account only for PCV, not for hib
delta_pne_pcv	<-	dat$pneumonia * pcv_effcov * data.predict$ppnedeaths_pcv
pne_pcv			<-	dat$pneumonia - delta_pne_pcv

delta_pne_hib	<- 	dat$pneumonia* data.predict$hib3/100 * data.predict$pne_hib_efficacy * data.predict$ppnedeaths_hib
pne_hib			<- 	dat$pneumonia- delta_pne_hib

delta_pne_pcv	<- 	dat$pneumonia* pcv_effcov * data.predict$ppnedeaths_pcv
pne_pcv			<- 	dat$pneumonia- delta_pne_pcv

delta_pne_hibpcv 	<- 	delta_pne_hib + delta_pne_pcv
pne_hibpcv 		<- 	dat$pneumonia- delta_pne_hibpcv

#* take account of the Hib and pcv coverage for meningitis
delta_men_hib	<- 	dat$meningitis* data.predict$hib3/100 * data.predict$men_hib_efsize * data.predict$pmendeaths_hib
men_hib			<- 	dat$meningitis- delta_men_hib

delta_men_pcv	<- 	dat$meningitis* pcv_effcov * data.predict$pmendeaths_pcv
men_pcv			<- 	dat$meningitis- delta_men_pcv

delta_men_hibpcv 	<- 	delta_men_hib + delta_men_pcv
men_hibpcv 		<- 	dat$meningitis- delta_men_hibpcv



#/* combination: pne_pcv	men_pcv	dia_rota
# pne_hib 		A
# pne_pcv 		B
# men_pcv 		C
# men_hib 		D
# dia_rota		E
#*/

#    cat("CHICK\n")
#    print(names(dat))
#    print(names(data.predict))

##############################################################################
#*permutation 1: ABC ACB
#*e.g. pne_hib pne_pcv men_pcv men_hib dia_rota
#*prop to allocate pne deaths averted by hib & pcv, allocate among 7 causes except pne
 pneumonia11 <- pne_hibpcv
    for( v in causes[which(causes != "pneumonia")] ) {
        dem<- rowSums(dat[,causes[which(causes !="pneumonia")]])
        dem<- ifelse(dem==0,1,dem)
        assign(paste0("p_",v,"11"), dat[,v] / dem )


        assign(paste0(v,"11") , dat[,v] + delta_pne_hibpcv * get(paste0("p_",v,"11")) )
}

if(verbose) print("PH 11")

if(verbose) {
  print(summary(data.predict$envelope_aids_measles_free))
  print(summary(pneumonia11))
}

 temp<- data.predict$envelope_aids_measles_free - pneumonia11 - injuries11 - malaria11 - meningitis11 -other11 - diarrhea11 - neonatal11 - congenital11
if(verbose){
  print(summary(temp))
    print(dim(data.predict))
}
    xt<-which(is.na(temp))
    if(verbose){
    print(xt)
    print(data.predict[which(is.na(temp)),])
    print(pneumonia11[xt])
    print(injuries11[xt])
    print(malaria11[xt])
    print(meningitis11[xt])
    print(other11[xt])
    print(diarrhea11[xt])
    print(neonatal11[xt])
    print(congenital11[xt])
}

    if(max(abs(temp)) > 1) {
            cat("\ncheck 11- should be 0 \n\n")
            print(summary(temp) )
            }

#    cat("\nCHICK 1x\n")
#*prop to allocate men deaths averted by hib & pcv, allocate among 7 cods except men
#gen delta_men_hib12	=	men11 * hib/100 * men_hib_efsize * pmendeaths_hib
 delta_men_hib12	<-	meningitis11 * data.predict$hib3/100 * data.predict$men_hib_efsize * data.predict$pmendeaths_hib
 delta_men_pcv12	<-	meningitis11 * pcv_effcov * data.predict$pmendeaths_pcv
 delta_men_hibpcv12 	<-	delta_men_hib12 + delta_men_pcv12
meningitis12 <- meningitis11 - delta_men_hibpcv12

for(v in causes[which(causes != "meningitis")]) {
        dem<- (pneumonia11+injuries11+malaria11+other11+diarrhea11+neonatal11+congenital11)
        dem<- ifelse(dem==0,1,dem)

    assign(paste0("p_",v,"12"), get(paste0(v, "11")) / dem )
assign(paste0(v, "12"), get(paste0(v,"11")) + delta_men_hibpcv12 * get(paste0("p_",v,"12")))
}
temp<- data.predict$envelope_aids_measles_free - pneumonia12 - injuries12 - malaria12 - meningitis12 -other12 - diarrhea12 - neonatal12 - congenital12
if(verbose) print("PH 12")
#print(summary(temp))
#print(summary(data.predict$envelope_aids_measles_free))


        if(max(abs(temp)) > 1) {
            cat("\ncheck 12- should be 0 \n\n")
    print(summary(temp) )
}

#*prop to allocate dia deaths averted by rota, allocate among 7 causes except dia
delta_dia_rota13 <- diarrhea12 * (data.predict$rota/100) * data.predict$dia_rota_effective * data.predict$pdiadeaths_rota
diarrhea13  <- diarrhea12 - delta_dia_rota13
    for( v in causes[which(causes != "diarrhea")]) {
        dem<-(pneumonia12+injuries12+malaria12+other12+meningitis12+neonatal12+congenital12)
        dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"13"), get(paste0(v,"12")) / dem)
assign(paste0(v,"13"), get(paste0(v,"12")) + delta_dia_rota13 * get(paste0("p_",v,"13")) )
}
temp<- data.predict$envelope_aids_measles_free - pneumonia13 - injuries13 - malaria13 - meningitis13 -other13 - diarrhea13 - neonatal13 - congenital13
if(verbose) print("PH 13")
        if(max(abs(temp)) > 1) {
            cat("\ncheck 13- should be 0 \n\n")
    print(summary(temp) )
}


#*******************************************************************************************************
#*permutation 2: ABECD ABEDC BAECD BAEDC
#*e.g. pne_hib pne_pcv dia_rota men_pcv men_hib

#*prop to allocate pne deaths averted by hib & pcv, allocate among 7 causes except pne
pneumonia21 <- pne_hibpcv
    for( v in causes[which(causes !="pneumonia")]) {
                dem<- rowSums(dat[,causes[which(causes !="pneumonia")]])
                dem<- ifelse(dem==0,1,dem)

                assign(paste0("p_",v,"21"), dat[,v] / dem )
                assign(paste0(v,"21"), dat[,v] + delta_pne_hibpcv * get(paste0("p_",v,"21")) )
}
temp<- data.predict$envelope_aids_measles_free - pneumonia21 - injuries21 - malaria21 - meningitis21 -other21 - diarrhea21 - neonatal21 - congenital21
if(verbose) print("PH 21")
        if(max(abs(temp)) > 1) {
            cat("\ncheck 21- should be 0 \n\n")
            print(summary(temp) )
}
#*prop to allocate dia deaths averted by rota, allocate among 7 causes except dia
delta_dia_rota22 <- diarrhea21 * (data.predict$rota/100) * data.predict$dia_rota_effective * data.predict$pdiadeaths_rota
diarrhea22  <- diarrhea21 - delta_dia_rota22
    for( v in causes[which(causes != "diarrhea")]) {
        dem<- (pneumonia21+injuries21+malaria21+other21+meningitis21+neonatal21+congenital21)
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"22"), get(paste0(v,"21")) / dem )
assign(paste0(v,"22"), get(paste0(v,"21")) + delta_dia_rota22 * get(paste0("p_",v,"22")) )
}

temp<- data.predict$envelope_aids_measles_free - pneumonia22 - injuries22 - malaria22 - meningitis22 -other22 - diarrhea22 - neonatal22 - congenital22
if(verbose) print("PH 22")
        if(max(abs(temp)) > 1) {
            cat("\ncheck 22- should be 0 \n\n")
    print(summary(temp) )
}


                                        #*prop to allocate men deaths averted by hib & pcv, allocate among 7 cods except men
delta_men_hib23	<-	meningitis22 * data.predict$hib3/100 * data.predict$men_hib_efsize * data.predict$pmendeaths_hib
delta_men_pcv23	<-	meningitis22 * pcv_effcov * data.predict$pmendeaths_pcv
delta_men_hibpcv23 	<- delta_men_hib23 + delta_men_pcv23
    meningitis23 <- meningitis22 - delta_men_hibpcv23
    if(verbose) print("PH 23")
    for( v in causes[which(causes != "meningitis")]) {
        dem<- (pneumonia22+injuries22+malaria22+other22+diarrhea22+neonatal22+congenital22)
        dem<- ifelse(dem==0,1,dem)
        assign(paste0("p_",v,"23"), get(paste0(v,22)) / dem)
        assign(paste0(v,"23"), get(paste0(v,22)) + delta_men_hibpcv23 * get(paste0("p_",v,23)) )
}
temp<- data.predict$envelope_aids_measles_free - pneumonia23 - injuries23 - malaria23 - meningitis23 -other23 - diarrhea23 - neonatal23 - congenital23
        if(max(abs(temp)) > 1) {
            cat("\ncheck 23- should be 0 \n\n")
    print(summary(temp) )
}
#*******************************************************************************************************
#*permutation 3: CDABE CDBAE DCABE DCBAE
#*e.g. men_pcv men_hib pne_hib pne_pcv dia_rota

#*prop to allocate men deaths averted by hib & pcv, allocate among 7 cods except men
#delta_men_pcv31	<-	dat$meningitis * pcv_effcov * data.predict$pmendeaths_pcv
                                        #meningitis31 <- dat$meningitis - delta_men_pcv31
#cat("\nmen_pcv\n")
#    print(summary(men_pcv))
delta_men_hib31	<-	dat$meningitis * data.predict$hib3/100 * data.predict$men_hib_efsize * data.predict$pmendeaths_hib
delta_men_pcv31	<-	dat$meningitis * pcv_effcov * data.predict$pmendeaths_pcv
delta_men_hibpcv31 <- delta_men_hib31 + delta_men_pcv31
    meningitis31<- dat$meningitis - delta_men_hibpcv31
    if(verbose) print("PH 31")
for( v in causes[which(causes != "meningitis")]) { # local vlpermmen {
        dem<-rowSums(dat[, causes[which(causes !="meningitis")] ])
        dem<- ifelse(dem==0,1,dem)
        assign(paste0("p_",v,"31"), dat[,v] / dem )
        assign(paste0(v,31), dat[,v] + delta_men_hibpcv31 * get(paste0("p_",v,31)) )
}
temp<- data.predict$envelope_aids_measles_free - pneumonia31 - injuries31 - malaria31 - meningitis31 -other31 - diarrhea31 - neonatal31 - congenital31
        if(max(abs(temp)) > 1) {
            cat("\ncheck 31- should be 0 \n\n")
            print(summary(temp) )
            }

#cat("\n\n CHICK 2\n")

#*prop to allocate pne deaths averted by hib & pcv, allocate among 7 causes except pne
delta_pne_hib32	<-	pneumonia31 * data.predict$hib3/100 * data.predict$pne_hib_efficacy * data.predict$ppnedeaths_hib
delta_pne_pcv32	<-	pneumonia31 * pcv_effcov * data.predict$ppnedeaths_pcv
delta_pne_hibpcv32 	<- delta_pne_hib32 + delta_pne_pcv32
    pneumonia32 <- pneumonia31 - delta_pne_hibpcv32
    if(verbose) print("PH 32")
    for( v in causes[which(causes != "pneumonia")]) {
            dem<-(injuries31+malaria31+meningitis31+other31+diarrhea31+neonatal31+congenital31)
            dem<- ifelse(dem==0,1,dem)
            assign(paste0("p_",v,"32"), get(paste0(v,31)) / dem)
            assign(paste0(v,32), get(paste0(v,"31")) + delta_pne_hibpcv32 * get(paste0("p_",v,32)) )
}
temp<- data.predict$envelope_aids_measles_free - pneumonia32 - injuries32 - malaria32 - meningitis32 -other32 - diarrhea32 - neonatal32 - congenital32
        if(max(abs(temp)) > 1) {
            cat("\ncheck 32- should be 0 \n\n")
print(summary(temp) )
}

#*prop to allocate dia deaths averted by rota, allocate among 7 causes except dia
delta_dia_rota33 <- diarrhea32 * (data.predict$rota/100) * data.predict$dia_rota_effective * data.predict$pdiadeaths_rota
    diarrhea33  <- diarrhea32 - delta_dia_rota33
    if(verbose) print("PH 33")
    for( v in  causes[which(causes != "diarrhea")]) {
            dem<- (pneumonia32+injuries32+malaria32+other32+meningitis32+neonatal32+congenital32)
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"33"), get(paste0(v,"32")) / dem )
assign(paste0(v,33), get(paste0(v,"32")) + delta_dia_rota33 * get(paste0("p_",v,33)) )
}
temp<- data.predict$envelope_aids_measles_free - pneumonia33 - injuries33 - malaria33 - meningitis33 -other33 - diarrhea33 - neonatal33 - congenital33
        if(max(abs(temp)) > 1) {
            cat("\ncheck 33- should be 0 \n\n")
    print(summary(temp) )
}
#*******************************************************************************************************
#*permutation 4: CDEAB CDEBA DCEAB DCEBA
#*e.g. men_pcv men_hib dia_rota pne_hib pne_pcv

#*prop to allocate men deaths averted by hib & pcv, allocate among 7 cods except men
#delta_men_pcv41	<-	dat$meningitis * pcv_effcov * data.predict$pmendeaths_pcv
delta_men_hib41	<-	dat$meningitis * data.predict$hib3/100 * data.predict$men_hib_efsize * data.predict$pmendeaths_hib
delta_men_pcv41	<-	dat$meningitis * pcv_effcov * data.predict$pmendeaths_pcv
delta_men_hibpcv41 	<-	delta_men_hib41 +delta_men_pcv41
    meningitis41 <- dat$meningitis - delta_men_hibpcv41
    if(verbose) print("PH 41")
#meningitis41 <- men_pcv #dat$meningitis - delta_men_pcv41
    for( v in causes[which(causes != "meningitis")]) {
            dem<- rowSums(dat[,causes[which(causes !="meningitis")]])
            dem<- ifelse(dem==0,1,dem)
            assign(paste0("p_",v,"41"), dat[,v] / dem )
            assign(paste0(v,41), dat[,v] + delta_men_hibpcv41 * get(paste0("p_",v,41)) )
}
#cat("\ncheck 41 - should be 0 \n\n")
temp<- data.predict$envelope_aids_measles_free - pneumonia41 - injuries41 - malaria41 - meningitis41 -other41 - diarrhea41 - neonatal41 - congenital41
        if(max(abs(temp)) > 1) {
            cat("\ncheck 41- should be 0 \n\n")
print(summary(temp) )
}
#*prop to allocate dia deaths averted by rota, allocate among 7 causes except dia
delta_dia_rota42 <- diarrhea41 * (data.predict$rota/100) * data.predict$dia_rota_effective * data.predict$pdiadeaths_rota
    diarrhea42  <- diarrhea41 - delta_dia_rota42
        if(verbose) print("PH 42")
    for( v in  causes[which(causes !="diarrhea")]) {
        dem<- (pneumonia41+injuries41+malaria41+other41+meningitis41+neonatal41+congenital41)
            dem<- ifelse(dem==0,1,dem)
        assign(paste0("p_",v,"42"), get(paste0(v,41)) / dem )
        assign(paste0(v,42), get(paste0(v,41)) + delta_dia_rota42 * get(paste0("p_",v,42)) )
}

#cat("\ncheck 42 - should be 0 \n\n")
    temp<- data.predict$envelope_aids_measles_free - pneumonia42 - injuries42 - malaria42 - meningitis42 -other42 - diarrhea42 - neonatal42 - congenital42
        if(max(abs(temp)) > 1) {
            cat("\ncheck 42- should be 0 \n\n")
    print(summary(temp) )
}
#*prop to allocate pne deaths averted by hib & pcv, allocate among 7 causes except pne
    delta_pne_hib43	<-	pneumonia42 * data.predict$hib3/100 * data.predict$pne_hib_efficacy * data.predict$ppnedeaths_hib
    delta_pne_pcv43	<-	pneumonia42 * pcv_effcov * data.predict$ppnedeaths_pcv
    delta_pne_hibpcv43 	<- delta_pne_hib43 + delta_pne_pcv43
    pneumonia43 <- pneumonia42 - delta_pne_hibpcv43
    if(verbose) print("PH 43")
    for( v in causes[which(causes !="pneumonia")]) {
        dem<- (injuries42+malaria42+meningitis42+other42+diarrhea42+neonatal42+congenital42)
            dem<- ifelse(dem==0,1,dem)
        assign(paste0("p_",v,"43"), get(paste0(v,42)) / dem )
        assign(paste0(v,43), get(paste0(v,42)) + delta_pne_hibpcv43 * get(paste0("p_",v,43)) )
}
#cat("\ncheck 43- should be 0 \n\n")
    temp<- data.predict$envelope_aids_measles_free - pneumonia43 - injuries43 - malaria43 - meningitis43 -other43 - diarrhea43 - neonatal43 - congenital43
        if(max(abs(temp)) > 1) {
            cat("\ncheck 43- should be 0 \n\n")
    print(summary(temp) )
        }


#*******************************************************************************************************
#*permutation 5: EABCD EABDC EBACD EBADC
#*e.g. dia_rota pne_hib pne_pcv men_pcv men_hib

#*prop to allocate dia deaths averted by rota, allocate among 7 causes except dia
delta_dia_rota51 <- dat$diarrhea * (data.predict$rota/100) * data.predict$dia_rota_effective * data.predict$pdiadeaths_rota
    diarrhea51  <- dat$diarrhea - delta_dia_rota51
        if(verbose) print("PH 51")
    for( v in  causes[which(causes != "diarrhea")]) {
        dem<- rowSums(dat[,causes[which(causes !="diarrhea")]])
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"51"), dat[,v] / dem )
assign(paste0(v,51), dat[,v] + delta_dia_rota51 * get(paste0("p_",v,51)) )
}
#cat("\ncheck 51- should be 0 \n\n")
    temp<- data.predict$envelope_aids_measles_free - pneumonia51 - injuries51 - malaria51 - meningitis51 -other51 - diarrhea51 - neonatal51 - congenital51
        if(max(abs(temp)) > 1) {
            cat("\ncheck 51- should be 0 \n\n")
    print(summary(temp) )
}
#*prop to allocate pne deaths averted by hib & pcv, allocate among 7 causes except pne
    delta_pne_hib52	<-	pneumonia51 * data.predict$hib3/100 * data.predict$pne_hib_efficacy * data.predict$ppnedeaths_hib
    delta_pne_pcv52	<-	pneumonia51 * pcv_effcov * data.predict$ppnedeaths_pcv
    delta_pne_hibpcv52 	<- delta_pne_hib52 + delta_pne_pcv52
    pneumonia52 <- pneumonia51 - delta_pne_hibpcv52


    for( v in causes[which(causes != "pneumonia")]) {
        dem<- (injuries51+malaria51+meningitis51+other51+diarrhea51+neonatal51+congenital51)
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"52"), get(paste0(v,51)) / dem )
assign(paste0(v,52), get(paste0(v,51)) + delta_pne_hibpcv52 * get(paste0("p_",v,52)) )
}
#cat("\ncheck 52 - should be 0 \n\n")
    temp<- data.predict$envelope_aids_measles_free - pneumonia52 - injuries52 - malaria52 - meningitis52 -other52 - diarrhea52 - neonatal52 - congenital52
        if(max(abs(temp)) > 1) {
            cat("\ncheck 52- should be 0 \n\n")
    print(summary(temp) )
}

#*prop to allocate men deaths averted by hib & pcv, allocate among 7 cods except men
#gen delta_men_hib53	=	men52 * hib/100 * men_hib_efsize * pmendeaths_hib
    delta_men_hib53	<-	meningitis52 * data.predict$hib3/100 * data.predict$men_hib_efsize * data.predict$pmendeaths_hib
    delta_men_pcv53	<-	meningitis52 * pcv_effcov * data.predict$pmendeaths_pcv
    delta_men_hibpcv53 	<- delta_men_hib53 + delta_men_pcv53
    meningitis53 <- meningitis52 - delta_men_hibpcv53

    for( v in  causes[which(causes !="meningitis")]) {
        dem<- (pneumonia52+injuries52+malaria52+other52+diarrhea52+neonatal52+congenital52)
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"53"), get(paste0(v,52)) / dem )
assign(paste0(v,53), get(paste0(v,52)) + delta_men_hibpcv53 * get(paste0("p_",v,53)) )
}
#cat("\ncheck 53 - should be 0 \n\n")
    temp<- data.predict$envelope_aids_measles_free - pneumonia53 - injuries53 - malaria53 - meningitis53 -other53 - diarrhea53 - neonatal53 - congenital53
            if(max(abs(temp)) > 1) {
            cat("\ncheck 53- should be 0 \n\n")
    print(summary(temp) )
}

#*******************************************************************************************************
#*permutation 6: ECDAB ECDBA EDCAB EDCBA
#*e.g. dia_rota men_pcv men_hib pne_hib pne_pcv

#*prop to allocate dia deaths averted by rota, allocate among 7 causes except dia
delta_dia_rota61 <- dat$diarrhea * (data.predict$rota/100) * data.predict$dia_rota_effective * data.predict$pdiadeaths_rota
diarrhea61  <- dat$diarrhea - delta_dia_rota61
    for( v in  causes[which(causes !="diarrhea")]) {
        dem<-rowSums(dat[,causes[which(causes !="diarrhea")]])
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"61"), dat[,v] / dem )
assign(paste0(v,61), dat[,v]  + delta_dia_rota61 * get(paste0("p_",v,61)) )
}
#cat("\ncheck 61- should be 0 \n\n")
    temp<- data.predict$envelope_aids_measles_free - pneumonia61 - injuries61 - malaria61 - meningitis61 -other61 - diarrhea61 - neonatal61 - congenital61
        if(max(abs(temp)) > 1) {
            cat("\ncheck 61- should be 0 \n\n")
    print(summary(temp) )
}

#*prop to allocate men deaths averted by hib & pcv, allocate among 7 cods except men
#gen delta_men_hib62	=	men61 * hib/100 * men_hib_efsize * pmendeaths_hib
    delta_men_hib62	<-	meningitis61 * data.predict$hib3/100 * data.predict$men_hib_efsize * data.predict$pmendeaths_hib
    delta_men_pcv62	<-	meningitis61 * pcv_effcov * data.predict$pmendeaths_pcv
    delta_men_hibpcv62 	<- delta_men_hib62 + delta_men_pcv62
    meningitis62 <- meningitis61 - delta_men_hibpcv62


#    delta_men_pcv62	<-	meningitis61 * pcv_effcov * data.predict$pmendeaths_pcv
#egen delta_men_hibpcv62 	=	rowtotal(delta_men_hib62 delta_men_pcv62)
#meningitis62 <- meningitis61 - delta_men_pcv62
    for( v in  causes[which(causes !="meningitis")]) {
        dem<- (pneumonia61+injuries61+malaria61+other61+diarrhea61+neonatal61+congenital61)
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"62"), get(paste0(v,61)) / dem )
assign(paste0(v,62), get(paste0(v,61)) + delta_men_hibpcv62 * get(paste0("p_",v,62)) )
}

temp<- data.predict$envelope_aids_measles_free - pneumonia62 - injuries62 - malaria62 - meningitis62 -other62 - diarrhea62 - neonatal62 - congenital62
        if(max(abs(temp)) > 1) {
#            cat("\ncheck 63- should be 0 \n\n")
    cat("\ncheck 62 - should be 0 \n\n")
    print(summary(temp) )
}


#*prop to allocate pne deaths averted by hib & pcv, allocate among 7 causes except pne
#gen delta_pne_hib63	=	pne62 * hib/100 * pne_hib_efficacy * ppnedeaths_hib
delta_pne_hib63	<-	pneumonia62 * data.predict$hib3/100 * data.predict$pne_hib_efficacy * data.predict$ppnedeaths_hib
delta_pne_pcv63	<-	pneumonia62 * pcv_effcov * data.predict$ppnedeaths_pcv
delta_pne_hibpcv63 	<- delta_pne_hib63 +delta_pne_pcv63
pneumonia63 <- pneumonia62 - delta_pne_hibpcv63


    for( v in  causes[which(causes !="pneumonia")]) {
            dem<- (injuries62+malaria62+meningitis62+other62+diarrhea62+neonatal62+congenital62)
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"63"), get(paste0(v,62)) / dem)
assign(paste0(v,63), get(paste0(v,62)) + delta_pne_hibpcv63 * get(paste0("p_",v,63)) )
}

    temp<- data.predict$envelope_aids_measles_free - pneumonia63 - injuries63 - malaria63 - meningitis63 -other63 - diarrhea63 - neonatal63 - congenital63
        if(max(abs(temp)) > 1) {
            cat("\ncheck 63- should be 0 \n\n")
            print(summary(temp) )
            }


#********************************************************************************************
#*obtain the average of the 6 permutations
#local vlperm0 pne inj mal men oth dia neo con
for( v in causes) {
assign(paste0(v,"_new"), (get(paste0(v,13))+get(paste0(v,23))+get(paste0(v,33))+
                            get(paste0(v,43))+get(paste0(v,53))+get(paste0(v,63)) )/6 )
}

    for( v in causes) {
        dem<-(pneumonia_new + injuries_new + malaria_new +
                                                         meningitis_new + other_new + diarrhea_new + neonatal_new + congenital_new)
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p",v,"_new"), get(paste0(v,"_new")) / dem)
assign(paste0(v,"_new"), get(paste0("p",v,"_new")) * data.predict$envelope_aids_measles_free )
}

    dem<- data.predict$envelope_aids_measles_free
    dem<- ifelse(dem==0,1,dem)
    pnorm <- (pneumonia_new + injuries_new + malaria_new +
           meningitis_new + other_new + diarrhea_new + neonatal_new + congenital_new) / dem
#print(summary(pnorm))



    temp<- data.predict$envelope_aids_measles_free - pneumonia_new - injuries_new - malaria_new - meningitis_new -other_new - diarrhea_new - neonatal_new - congenital_new
    if(max(abs(temp)) > 1) {
    cat("\ncheck new - should be 0 \n\n")
    print(summary(temp) )
    }

#rename (ppne_new pinj_new pmal_new pmen_new poth_new pdia_new pneo_new pcon_new) (p0_new p1_new p2_new p3_new p4_new p5_new p6_new p7_new)

#********************************************************************************************
# *assign 0 to countries where pfalciparum is not present -- pfpr==0
#    malaria_new_redis_final<- malaria_new
# replace mal_new_redis_final=0 if pfpr==0
    malaria_new_redis_final<- ifelse(data.predict$pfpr==0, 0, malaria_new)
    delta_mal<- malaria_new-malaria_new_redis_final

    for( v in  causes[which(causes != "malaria")]) {
        dem<-(pneumonia_new+injuries_new+meningitis_new+other_new+diarrhea_new+neonatal_new+congenital_new)
            dem<- ifelse(dem==0,1,dem)
assign(paste0("p_",v,"_mal"), get(paste0(v,"_new")) / dem )
assign(paste0(v,"_new_redis_final"), get(paste0(v,"_new")) + get(paste0("p_",v,"_mal") )*delta_mal)
}



    for( v in causes) {
            dem<- data.predict$envelope_aids_measles_free
            dem<- ifelse(dem==0,1,dem)
            assign(paste0("p_",v),get(paste0(v,"_new_redis_final")) / dem)
 }



dat$pneumonia<-p_pneumonia
dat$malaria<- p_malaria
dat$other<- p_other
dat$meningitis<- p_meningitis
dat$injuries<- p_injuries
dat$diarrhea<- p_diarrhea
dat$neonatal<- p_neonatal
dat$congenital<- p_congenital

    dat<- dat[,causes]
    if( max(abs( rowSums(dat) - 1)) >0.0001) {
    cat("\ncheck - should be 1 \n\n")
    print(summary(rowSums(dat)))
    print(rownames(dat[which(rowSums(dat)==0),]))
    }
    return(dat)
}

