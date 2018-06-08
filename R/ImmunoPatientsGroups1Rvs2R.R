#################################
# this section is to define groups of patients treated with ICB that are primary resistant or secondary resitant and the time of biospy with respect to treatment
# 1R is primary resistant
# 2R is secondary resistant

# PFS discretisation: define the cutoff for resistant, based on the litterature
ImmunoMoscatoPatients$BadPFS <- ImmunoMoscatoPatients$PFS_immuno<90
ImmunoMoscatoPatients$GoodPFS <- ImmunoMoscatoPatients$PFS_immuno>90

# groups
NUM2Rpost<- ImmunoMoscatoPatients$Personid[ImmunoMoscatoPatients$BiopsyAfterProg&ImmunoMoscatoPatients$GoodPFS]
NUM2Rpost<-NUM2Rpost[!duplicated(NUM2Rpost)]
NUM2Rpost<-NUM2Rpost[!is.na(NUM2Rpost)]

NUM1Rpost<- ImmunoMoscatoPatients$Personid[ImmunoMoscatoPatients$BiopsyAfterProg&ImmunoMoscatoPatients$BadPFS]
NUM1Rpost<-NUM1Rpost[!duplicated(NUM1Rpost)]
NUM1Rpost<-NUM1Rpost[!is.na(NUM1Rpost)]
NUM1Rpost
length(NUM1Rpost)

NUM2Rpre<- ImmunoMoscatoPatients$Personid[ImmunoMoscatoPatients$BiopsyBeforeInclusion&ImmunoMoscatoPatients$GoodPFS&ImmunoMoscatoPatients$RNAdata]
NUM2Rpre<-NUM2Rpre[!duplicated(NUM2Rpre)&!is.na(NUM2Rpre)]

NUM1Rpre<- ImmunoMoscatoPatients$Personid[ImmunoMoscatoPatients$BiopsyBeforeInclusion&ImmunoMoscatoPatients$BadPFS&ImmunoMoscatoPatients$RNAdata]
NUM1Rpre<-NUM1Rpre[!duplicated(NUM1Rpre)]
NUM1Rpre<-NUM1Rpre[!is.na(NUM1Rpre)]

DUPL<-c(NUM1Rpost,NUM2Rpost)[c(NUM1Rpost,NUM2Rpost)%in%c(NUM1Rpre , NUM2Rpre)]

# another way to constitute groups are:
BEFORE<-ImmunoMoscatoPatients$BiopsyBeforeInclusion&ImmunoMoscatoPatients$RNAdata
NUMRpre<-ImmunoMoscatoPatients$Personid[BEFORE]
NUMRpre<-NUMRpre[!is.na(NUMRpre)]
NUMRpre<-NUMRpre[!duplicated(NUMRpre)]

NOTBEFORE<-!ImmunoMoscatoPatients$BiopsyBeforeInclusion&ImmunoMoscatoPatients$RNAdata
NUMRpost<-ImmunoMoscatoPatients$Personid[NOTBEFORE]
NUMRpost<-NUMRpost[!is.na(NUMRpost)]
NUMRpost<-NUMRpost[!duplicated(NUMRpost)]

# check patients that are in both groups (that have either a biopsy before and after or in-between a treatment
DUPL<-NUMRpre[NUMRpre%in%NUMRpost]

#keep these groups in memory, you'll need it furhter on
