##################
# Data formats processing
# author: Loic Verlingue, DITEP
##################

##################
# 1. Gene expression computed with Salmon (first version on November 3rd 2017)
# txtimport data matrices stored in "SalmonGENEXP"

## for counts
# clean samplle id
NUMFILE<-gsub("_.*","", gsub("foie","", gsub("RE","", gsub("PED","", gsub("Sample_","", gsub("-.*", "", colnames(SalmonGENEXP$counts)))))))
colnames(SalmonGENEXP$counts)<-NUMFILE

#turn into log10 transformed counts
SalmonCounts<-log10(SalmonGENEXP$counts+1)

## for tpm, idem
NUMFILE<-gsub("_.*","", gsub("foie","", gsub("RE","", gsub("PED","", gsub("Sample_","",gsub("-.*", "", colnames(SalmonGENEXP$abundance)))))))
colnames(SalmonGENEXP$abundance)<-NUMFILE
SalmonTPM<-log10(SalmonGENEXP$abundance+1)

##################
# 2. clinical data, from patients' table called "ImmunoMoscatoPatients" (first version on November 14th 2017)
# clean sample id and tag those with RNAseq
NUMFILE<-gsub("_.*","", gsub("foie","", gsub("RE","", gsub("PED","", gsub("Sample_","", gsub("-.*", "", colnames(SalmonGENEXP$counts)))))))
ImmunoMoscatoPatients$RNAdata <- ImmunoMoscatoPatients$Personid%in%NUMFILE

# tag time of biopsy relative to time of C1D1 of treatment
ImmunoMoscatoPatients$BiopsyBeforeInclusion<-ImmunoMoscatoPatients$date_biopsie-ImmunoMoscatoPatients$date_C1D1_immuno<0
ImmunoMoscatoPatients$BiopsyAfterProg <- ImmunoMoscatoPatients$date_biopsie-ImmunoMoscatoPatients$date_progression_Recist>-20
ImmunoMoscatoPatients$Interval_PD_biops<-as.numeric(ImmunoMoscatoPatients$date_biopsie-ImmunoMoscatoPatients$date_progression_Recist)

# compute PFS, OS using the Surv function from survival package
# and compute ORR as:
ImmunoMoscatoPatients$ORR[ImmunoMoscatoPatients$RR_ML%in%c("PD","SD")]<-"PD_SD"
ImmunoMoscatoPatients$ORR[ImmunoMoscatoPatients$RR_ML%in%c("PR","CR")]<-"PR_CR"

# limit skewed variables for tumor types and site of biopsy
ImmunoMoscatoPatients[ImmunoMoscatoPatients$GlobalHisto%in%c("","CONNECTIVE AND SOFT TISSUS","GLANDS","SNC","SKIN"),"GlobalHisto"]<-"OTHERS"
ImmunoMoscatoPatients[ImmunoMoscatoPatients$Loc_biopsie%in%c("","biop derm","mediastin","SNC","SKIN"),"Loc_biopsie"]<-"autre"
ImmunoMoscatoPatients[ImmunoMoscatoPatients$Loc_biopsie%in%"biop derm","Loc_biopsie"]<-"peau"
ImmunoMoscatoPatients[ImmunoMoscatoPatients$Loc_biopsie%in%c("peritoine","surrenale", "pancreas","vesicale","rate","rein","retro-peritoneal"),"Loc_biopsie"]<-"abdominal"
ImmunoMoscatoPatients[ImmunoMoscatoPatients$Loc_biopsie%in%c("prostate","vagin"),"Loc_biopsie"]<-"pelvis"
ImmunoMoscatoPatients$Loc_biopsie<-gsub(" ","", ImmunoMoscatoPatients$Loc_biopsie)

# further on, a step of cleaning biopsy variable
ImmunoMoscatoPatients$Loc_biopsie[ImmunoMoscatoPatients$Loc_biopsie==""|is.na(ImmunoMoscatoPatients$Loc_biopsie)]<-"others"
ImmunoMoscatoPatients$Loc_biopsie[ImmunoMoscatoPatients$Loc_biopsie=="poumon"]<-"lung"
ImmunoMoscatoPatients$Loc_biopsie[ImmunoMoscatoPatients$Loc_biopsie=="foie"]<-"liver"
ImmunoMoscatoPatients$Loc_biopsie[ImmunoMoscatoPatients$Loc_biopsie=="ganglions"]<-"Lymph node"
ImmunoMoscatoPatients$Loc_biopsie[ImmunoMoscatoPatients$Loc_biopsie=="abdominal"]<-"abdomen"
ImmunoMoscatoPatients$Loc_biopsie[ImmunoMoscatoPatients$Loc_biopsie%in%c("orl","orl ")]<-"head & neck"


####################################
# 3.FusTable: a table that stores several usefull clinical information 
 
FusTable$RNAdata<-FusTable$Personid%in%NUMFILE

# discard duplicates
table(duplicated(FusTable$Personid))
FusTable<-FusTable[!duplicated(FusTable$Personid),]

# keep only patients with RNA data
table(FusTable$RNAdata)
FusTable<-FusTable[FusTable$RNAdata,]

# allocate miscellaneous or unknown
FusTable$TCGA[is.na(FusTable$TCGA)]<-"Unknown"
FusTable$TCGA[FusTable$TCGA==""]<-"miscellaneous"

###########################################
# 4. Quality: table of the quality check for patients' samples before sequencing

Quality$RIN_ARN<-as.numeric(Quality$RIN)
head(Quality$RIN);head(Quality$RIN_ARN)
Quality$NIP<-paste(substring(Quality$NumIGR,first = 0, last = 4),"-",substring(Quality$NumIGR, first = 5),sep = "")
Quality$Personid<-gsub(" .*","", gsub("foie","", gsub("RE","", gsub("PED","", gsub("_.*", "", Quality$Num.Incl)))))
Quality$NumRNA<- gsub(" ","", gsub("-ARN", "", Quality$Nom.Ech.1))

# date biopsy
Quality$date_biopsie<-as.Date(Quality$date.pvt.2, format="%d/%m/%Y")
Quality$date_biopsie[is.na(Quality$date_biopsie)]<-as.Date(Quality$date.pvt.1, format="%d/%m/%Y")[is.na(Quality$date_biopsie)]
Quality$date_biopsie[is.na(Quality$date_biopsie)]<-as.Date(Quality$date.pvt, format="%d/%m/%Y")[is.na(Quality$date_biopsie)]

# add date biopsy from fustable when not in Quality table
DATE_RNA<-FusTable[FusTable$Personid%in%Quality$Personid[is.na(Quality$date_biopsie)],c("Personid","date_biopsie")]
for(i in seq(nrow(DATE_RNA))){
  Quality[Quality$Personid%in%DATE_RNA$Personid[i]&is.na(Quality$date_biopsie),"date_biopsie"]<-DATE_RNA[i,2]
}

