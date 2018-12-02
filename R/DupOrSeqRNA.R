##############
# this script allows to identify samples that are duplicates of a single biopsy, and keep those that are sequencial biopsies of a single patient
# remark: I have used dates of biopsies and sample ID and intersected it from tables: SalmonGENEXPset, Quality and Fustable
# see DataFormatProcessing.R for descriptions of these tables
# author: Loic Verlingue
##############

# retrieve patients that have duplicated RNAseq analysis
# can be : re-biopsy, 2nd analysis (for poor quality in first analysis or other reasons). For the latter, check expression quality.

DUPL<-NUMFILE[duplicated(NUMFILE)]
DUPL<-DUPL[!duplicated(DUPL)] 

# correlation of duplicated patients expression with mean genes' expression in the cohort
# usefull to filter on quality, if this is the reason for multiple analysis on the same biopsy
RMEANS<-rowMeans(SalmonGENEXPset)
DupCor<-sapply(DUPL,function(x){
  return(cor(SalmonGENEXPset[,colnames(SalmonGENEXPset)[NUMFILE%in%x]],RMEANS))
})
rm(RMEANS)
DupCor
NAMESdup<-unlist(lapply(DupCor,rownames))

### when a patient has 2 analyis, discard duplicates but keep sequencial analysis
# for that, Use Quality table
table(NAMESdup%in%Quality$NumRNA)

# remove samples that are not in Quality and update NAMESdup
table(colnames(SalmonGENEXPset)%in%NAMESdup[!NAMESdup%in%Quality$NumRNA])
SalmonGENEXPset<-SalmonGENEXPset[,!colnames(SalmonGENEXPset)%in%NAMESdup[!NAMESdup%in%Quality$NumRNA]]
dim(SalmonGENEXPset);table(NAMESdup%in%colnames(SalmonGENEXPset))
NAMESdup<-NAMESdup[NAMESdup%in%colnames(SalmonGENEXPset)]

# remove samples that have same date of biopsy (real duplicates)
DUPS<-apply( Quality[match(NAMESdup,Quality$NumRNA),c("Num.Incl","date_biopsie")],1,function(x)paste(x,collapse = ""))
DUPS
  # filter those by best quality (most correlated to exp mat)
ToDiscard<-lapply(DupCor,function(x){
  x<-as.data.frame(x)
  DISC<-rownames(x)[!rownames(x)%in%NAMESdup]
  x<-x[rownames(x)%in%NAMESdup,,drop=F]
  if( any(duplicated(DUPS[NAMESdup%in%rownames(x)])) ){
    # remove the one less correlated 
    DISC<-c(DISC, rownames(x)[!x[,1]%in%max(x[,1])] ) 
  } 
  return(DISC)
})

ToDiscard<-as.character(unlist(ToDiscard))
table(colnames(SalmonGENEXPset)%in%ToDiscard)

SalmonGENEXPset<-SalmonGENEXPset[,!colnames(SalmonGENEXPset)%in%ToDiscard]
NUMFILE<-gsub("_.*","", gsub("foie","", gsub("RE","", gsub("PED","", gsub("Sample_","", gsub("-.*", "", colnames(SalmonGENEXPset)))))))
dim(SalmonGENEXPset)
