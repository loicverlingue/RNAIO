##################
# Working on data format
##################

##################
# 1. Gene expression computed with Salmon, processed on November 3rd 2017
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
# 2. clinical data, from patients' table called "ImmunoMoscatoPatients", updated in November 14th 2017
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
