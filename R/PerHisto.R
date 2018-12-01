#################
# numbers per tumor types and relative gene-based expression
# author: Loic Verlingue, DITEP
################

#### for variable names, refer to DataFormatProcessing.R 

### patients numbers
TCGAtable<-data.frame()
for(i in NUMFILE){
  if(i%in%FusTable$Personid){
    TCGAtable[i,"Personid"]<-i
    TCGAtable[i,"TCGA"]<-FusTable$TCGA[FusTable$Personid%in%i]
    TCGAtable[i,"GlobalHisto"]<-FusTable$GlobalHisto[FusTable$Personid%in%i]
    TCGAtable[i,"fusion"]<-FusTable$fusion[FusTable$Personid%in%i]
  }
}

sort(table(TCGAtable$TCGA))

##### expression of genes on the full database

pdf(paste(dir,"Expression_by_histo_log.pdf",sep = ""))

for(Gene in GENES){
  ORD1<-data.frame(row.names = unique(TCGAtable$TCGA))
  for(i in unique(TCGAtable$TCGA)){
    SET<-SalmonGENEXPset[rownames(SalmonGENEXPset)%in%Gene,NUMFILE%in%TCGAtable$Personid[TCGAtable$TCGA%in%i]]
    ORD1[i,"N"]<-length(SET)
    ORD1[i,"median"]<-median(SET,na.rm = T)
    ORD1[i,"mean"]<-mean(SET,na.rm = T)
    ORD1[i,"sd"]<-sd(SET,na.rm = T)
  }
  
  ORD1<-ORD1[order(ORD1$median,decreasing = T),]
  ORD1<-ORD1[c(rownames(ORD1[!is.na(ORD1$sd),]),rownames(ORD1[is.na(ORD1$sd),])),]

  COL<-rainbow(length(unique(TCGAtable$TCGA)))
  par(mar= c(8, 4, 4, 2) + 0.1)
  plot(1,1,type='n', axes=F, xlab="",xlim=c(0,length(unique(TCGAtable$TCGA))+1), ylim=c(0,12),
       main=paste(Gene,"expression"), ylab="log2(TPM+1)")
  y=0
  for(i in rownames(ORD1)){
    y=y+1
    SET<-SalmonGENEXPset[rownames(SalmonGENEXPset)%in%Gene,NUMFILE%in%TCGAtable$Personid[TCGAtable$TCGA%in%i]]
    boxplot(SET, at=y,col=COL[y],axes=F,add=T)  
  }
  
  axis(2,las=2)
  mtext(text = rownames(ORD1),side = 1,at = seq(y),las=2)
 write.csv2(ORD1,paste(dir, Gene, "Expression_by_histo.csv", sep = ""))
}
dev.off()
