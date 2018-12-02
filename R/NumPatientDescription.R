#####################
# author: Loic Verlingue, DITEP, 08/06/2018
# description of the number of patients with RNAseq and the type of treatment they have received (ICB or MTA) and respective time of biopsy
# this script implies you have already :
# - loaded RNAseq data (detailed in txtimport.R) 
# - processed with DataFromat.R to obtain the ImmunoMoscatoPatients table for patients treated with ICB 
# - done the groups of ICB treatment responses with ImmunoPatientsGroups1Rvs2R.R
#
# This section will load other tables: 
# Clinic - is the CRF of the MOSCATO-01 trial for patients treated with MTA
# Fustable - is a table with minimal characteristics (age, tumor type, time of biopsy...) for real time monitoring
# dir - is a directory I use
#####################

#load supplementary tables
# clinic is the moscato table - treatment with MTA
Clinic$NIP <- as.character(Clinic$NIP)
Clinic$OS <- as.Date(Clinic$dfu, format = "%d/%m/%Y")-as.Date(Clinic$icdtc, format = "%d/%m/%Y")
ORTR <- Clinic[Clinic$oriente_traite==1,]

# Processing Fustable
# tag those with RNAseq
FusTable$RNAdata<-FusTable$Personid%in%NUMFILE

# keep only patients with RNA data
table(FusTable$RNAdata)
FusTable<-FusTable[FusTable$RNAdata,]

# discard duplicates
FusTable<-FusTable[!duplicated(FusTable$Personid),]

# allocate miscellaneous or unknown
FusTable$TCGA[is.na(FusTable$TCGA)]<-"Unknown"
FusTable$TCGA[FusTable$TCGA==""]<-"miscellaneous"
FusTable$GlobalHisto[is.na(FusTable$GlobalHisto)]<-"Unknown"
FusTable$GlobalHisto[FusTable$GlobalHisto==""]<-"miscellaneous"

# to allocate treatment types, needs to define 1st the categories
#

# then integrate that into Fustable
FusTable[FusTable$Personid%in%Clinic$Personid[Clinic$oriente_traite==1],"TTT"]<-"MTA"
FusTable[FusTable$Personid%in%NUMRpost,"TTT"]<-"IONR"
FusTable[FusTable$Personid%in%NUMRpre,"TTT"]<-"IO"
FusTable[FusTable$Personid%in%DUPL,"TTT"]<-"both_IO_IONR"
FusTable[FusTable$Personid%in%IOMTA,"TTT"]<-"IO_NR_MTA"
FusTable$TTT[is.na(FusTable$TTT)]<-"Others"
table(FusTable$TTT)

# do a new table to describe the patients' number per tumor types (Histo) and treatments

# this one is for the organ of origin of the tumor
Histo<-as.data.frame(sort(table(FusTable[FusTable$Personid%in%NUMFILE,"GlobalHisto"]),decreasing = T))
for(i in Histo$Var1){
  Histo[Histo$Var1%in%i,"IONR"]<-nrow(FusTable[FusTable$GlobalHisto%in%i&FusTable$Personid%in%NUMRpost,])
  Histo[Histo$Var1%in%i,"IO"]<-nrow(FusTable[FusTable$GlobalHisto%in%i&FusTable$Personid%in%NUMRpre,])
  Histo[Histo$Var1%in%i,"MTA"]<-nrow(FusTable[FusTable$GlobalHisto%in%i&FusTable$Personid%in%Clinic$Personid[Clinic$oriente_traite==1],])
  Histo[Histo$Var1%in%i,"IO_IONR"]<-nrow(FusTable[FusTable$GlobalHisto%in%i&FusTable$Personid%in%DUPL,])
  Histo[Histo$Var1%in%i,"IO_NR_MTA"]<-nrow(FusTable[FusTable$GlobalHisto%in%i&FusTable$Personid%in%IOMTA,])
}

sum(Histo$Freq)
write.csv2(Histo, paste(dir,"/DescriptionCohort_Organ_GlobalHisto.csv",sep = ""))

# this one is for the tumor type with TCGA like nomenclature (column "TCGA")
Histo<-as.data.frame(sort(table(FusTable[FusTable$Personid%in%NUMFILE,"TCGA"]),decreasing = T))
for(i in Histo$Var1){
  Histo[Histo$Var1%in%i,"IONR"]<-nrow(FusTable[FusTable$TCGA%in%i&FusTable$Personid%in%NUMRpost,])
  Histo[Histo$Var1%in%i,"IO"]<-nrow(FusTable[FusTable$TCGA%in%i&FusTable$Personid%in%NUMRpre,])
  Histo[Histo$Var1%in%i,"MTA"]<-nrow(FusTable[FusTable$TCGA%in%i&FusTable$Personid%in%Clinic$Personid[Clinic$oriente_traite==1],])
  Histo[Histo$Var1%in%i,"IO_IONR"]<-nrow(FusTable[FusTable$TCGA%in%i&FusTable$Personid%in%DUPL,])
  Histo[Histo$Var1%in%i,"IO_NR_MTA"]<-nrow(FusTable[FusTable$TCGA%in%i&FusTable$Personid%in%IOMTA,])
}

sum(Histo$Freq)
write.csv2(Histo, paste(dir,"/DescriptionCohort_TCGAnomenclature.csv",sep = ""))
