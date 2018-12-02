##################
# author : Loic Verlingue, 08/06/2018
# this section allows to identify patients with sequencial, meaning pre post ICB biopsies
# has to be imporved for details 
# it requires to load another table - Quality - with the date of RNAseq reltively to the histo identifier, to avoid problems!
#################

# load fustable again
FusTable$IOttt<-FusTable$NIP%in%ImmunoMoscatoPatients$NIP

# load data Quality
Quality<-read.csv2("/Bilan_.csv", stringsAsFactors = F, dec = ".")

Quality$RIN_ARN<-as.numeric(Quality$RIN)
Quality$NIP<-paste(substring(Quality$Num_IGR,first = 0, last = 4),"-",substring(Quality$Num_IGR, first = 5),sep = "")
Quality$Personid<-gsub(" .*","", gsub("foie","", gsub("RE","", gsub("PED","", gsub("_.*", "", Quality$Num_Incl)))))
Quality$NumRNA<- gsub("-ARN", "", Quality$Nom_Ech.1)

# there may be a way to increase the matches but....
DATE_RNA<-Quality[Quality$NumRNA%in%gsub("_quant", "",colnames(SalmonGENEXP$abundance)),c("NumRNA","date_pvt.2")]
DATE_RNA<-DATE_RNA[!duplicated(DATE_RNA$NumRNA),]

# do a table with sequencial biopsies (maybe has to be improved for patients with same personid)
DUP<-FusTable[FusTable$Personid%in%ImmunoMoscatoPatients$Personid|FusTable$MATCHR%in%ImmunoMoscatoPatients$Personid|FusTable$MOSCATO%in%ImmunoMoscatoPatients$Personid,]
DUP<-FusTable[FusTable$NIP%in%ImmunoMoscatoPatients$NIP,]
DUP<-DUP[!is.na(DUP$NIP),]

# put the date of RNAseq analysis
for(i in 1:nrow(DUP)){
  for(y in colnames(DATE_RNA)){
    VAR<-DATE_RNA[DATE_RNA$NumRNA%in%gsub("_quant","", colnames(SalmonGENEXPset)[NUMFILE%in%DUP$Personid[i]]),y]
    if(length(VAR)>0){
      DUP[i,y]<-DATE_RNA[DATE_RNA$NumRNA%in%gsub("_quant","", colnames(SalmonGENEXPset)[NUMFILE%in%DUP$Personid[i]]),y]
    }
  }
}

# check date biopsy and date RNAseq is the same
DUP$date_pvt.2<-as.Date(DUP$date_pvt.2,format = "%d/%m/%Y")

DUP<-DUP[duplicated(DUP$NIP,fromLast = F)|duplicated(DUP$NIP,fromLast = T),]

# put in order
DUP<-DUP[order(DUP$date_biopsie),]
DUP<-DUP[order(DUP$NIP),]

x=DUP$NIP[1]
# check dates of biopsies are not the same
DUPRNA<-sapply(unique(DUP$NIP), function(x) {
  c(length(unique(DUP[DUP$NIP%in%x,"date_biopsie"])),
  length(unique(DUP[DUP$NIP%in%x,"RNAdata"])))
})
DUPRNA
DUP<-DUP[DUP$NIP%in%colnames(DUPRNA)[DUPRNA[2,]==1],]

DUPRNA<-sapply(unique(DUP$NIP), function(x) {
  unique(DUP[DUP$NIP%in%x,"RNAdata"])==T
})
DUPRNA
DUP<-DUP[DUP$NIP%in%names(DUPRNA)[DUPRNA==T],]
DUP[,c("NIP", "Personid","RNAdata","IOttt","date_biopsie")]

# check time of biopsy
ImmunoMoscatoPatients$date_progression_Recist
DUPIO<-ImmunoMoscatoPatients[ImmunoMoscatoPatients$NIP%in%DUP$NIP,c("NIP","date_C1D1_immuno","date_progression_Recist")]

# discard patients with 2 biopsies at the same time regarding administration of IO
SEQBIOPS<-sapply(unique(DUP$NIP), function(x){
  length(unique(DUP[DUP$NIP%in%x,"date_biopsie"]<DUPIO[DUPIO$NIP%in%x,"date_C1D1_immuno"]))
})
SEQBIOPS

DUP<-DUP[DUP$NIP%in%names(which(SEQBIOPS==2)),]
DUP[,c("NIP", "Personid","RNAdata","IOttt","date_biopsie")]

##########################
# do the plots

# reload RNAseq
SalmonGENEXPset<-log2(SalmonGENEXP$abundance+1)
DUP[,c("NIP", "Personid","RNAdata","IOttt","date_biopsie","date_pvt.2")]

#loop

pdf(paste(dir,"SeqBiops.pdf",sep = ""))
j=0
for(i in 1:6){ # 6 is the number of patients with sequencial biopsies in this analysis
  j=j+1
  MAT<-cbind(SalmonGENEXPset[GENES,grep(DUP$Personid[j],colnames(SalmonGENEXPset))],SalmonGENEXPset[GENES,grep(DUP$Personid[j+1],colnames(SalmonGENEXPset))[1]])
  COL<-rainbow(nrow(MAT))
  matplot(t(MAT), type = "l",col = COL,xlim=c(1,3),ylim = c(min(MAT)-1,max(MAT)) ,axes=F,ylab = "log2(TPM+1)",xlab = "",
          lty=1, lwd=2, main=paste("Patient",i) )
  legend('topright',col=COL,legend = rownames(MAT), lty=1,lwd=2)
  axis(2,las=2)
  
  DUR1<-as.numeric(diff.Date(c(DUP$date_biopsie[j],ImmunoMoscatoPatients[ImmunoMoscatoPatients$Personid%in%DUP$Personid[j],"date_C1D1_immuno"])))
  POINT1<-1+as.numeric(diff.Date(c(DUP$date_biopsie[j],ImmunoMoscatoPatients[ImmunoMoscatoPatients$Personid%in%DUP$Personid[j],"date_C1D1_immuno"])))/
    as.numeric(diff.Date(c(DUP$date_biopsie[j],DUP$date_biopsie[j+1])))
  DUR2<-as.numeric(diff.Date(ImmunoMoscatoPatients[ImmunoMoscatoPatients$Personid%in%DUP$Personid[j],c("date_C1D1_immuno","date_progression_Recist")]))
  POINT2<-POINT1+as.numeric(diff.Date(ImmunoMoscatoPatients[ImmunoMoscatoPatients$Personid%in%DUP$Personid[j],c("date_C1D1_immuno","date_progression_Recist")]))/
    as.numeric(diff.Date(c(DUP$date_biopsie[j],DUP$date_biopsie[j+1])))
  
  arrows(y0 = min(MAT)-1,y1 = min(MAT)-1,x0 = POINT1,x1 = POINT2,length = 0.1)
  mtext(text = "Duration of IO",side = 1,line = 0.5,at = mean(c(POINT1,POINT2)))
  mtext(text = paste("day",c(DUR1,DUR2)), side = 1,line = 0,at = c(POINT1,POINT2) )
  
  COR<-cor(SalmonGENEXPset[,grep(DUP$Personid[j],colnames(SalmonGENEXPset))],SalmonGENEXPset[,grep(DUP$Personid[j+1],colnames(SalmonGENEXPset))[1]])
  legend("bottomright", legend = paste("R² whole 2 biopsies",round(COR,4)), bty = "n")
  j=j+1
}
dev.off()
