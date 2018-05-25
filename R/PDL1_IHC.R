# IHC from Moscato

#pooling the data from July 3rd 2017
IHCPDL1<-read.csv2("/IHC/FinalTable tableIOMOSC_IHCPDL1_20170703.csv")
IHCPDL1<-IHCPDL1[!is.na(IHCPDL1$CTUM0),]

IHCPDL1N2<-read.csv2("/IHC/IHC_PDL1_MOSC.csv")
IHCPDL1N2<-IHCPDL1N2[!is.na(IHCPDL1N2$CTUM0),]

COLN<-colnames(IHCPDL1)[colnames(IHCPDL1)%in%colnames(IHCPDL1N2)]

IHCPDL1<-rbind(IHCPDL1[,COLN],IHCPDL1N2[,COLN])

# definition of PDL1 groups
# JY Scoazec proposes: 3gps = 0, 1-50%, >50%
IHCPDL1[IHCPDL1$CTUM0==100,"GP"]<-0
IHCPDL1[IHCPDL1$CTUM2.>=50|IHCPDL1$CTUM3.>=50,"GP"]<-2
IHCPDL1[is.na(IHCPDL1$GP),"GP"]<-1
IHCPDL1$GP

# Integrate in clinical table
COLN<-c("Target","PFS_immuno","SurvPFS_immuno","OS","SurvOS")
COLN%in%colnames(ImmunoMoscatoPatients)

for(i in seq(nrow(ImmunoMoscatoPatients)) ){
  if(length(table(IHCPDL1$Num_Incl_short%in%ImmunoMoscatoPatients$Personid[i]))==2){
    ImmunoMoscatoPatients[i,"IHCPDL1"]<-IHCPDL1[IHCPDL1$Num_Incl_short%in%ImmunoMoscatoPatients$Personid[i],"GP"][1]
  } else {
    print(ImmunoMoscatoPatients$Personid[i])
  }
}

NUM<-ImmunoMoscatoPatients$RNAdata&ImmunoMoscatoPatients$BiopsyBeforeInclusion

IHCNA<-ImmunoMoscatoPatients$IHCPDL1
NUMleft<-ImmunoMoscatoPatients$Personid[is.na(IHCNA)&ImmunoMoscatoPatients$RNAdata&ImmunoMoscatoPatients$BiopsyBeforeInclusion]
NUMIHC<-!is.na(IHCNA)&ImmunoMoscatoPatients$RNAdata&ImmunoMoscatoPatients$BiopsyBeforeInclusion

library(survcomp)
#par(mar=c(8,2,4,2)+0.1)
VALUESPFS<-summary(coxph(ImmunoMoscatoPatients$SurvPFS_immuno[NUMIHC]~ImmunoMoscatoPatients$IHCPDL1[NUMIHC]))
km.coxph.plot(formula.s=Surv(ImmunoMoscatoPatients$PFS_immuno[NUMIHC]/30,ImmunoMoscatoPatients$PD_immuno[NUMIHC])~ImmunoMoscatoPatients$IHCPDL1[NUMIHC], 
              mark.time=TRUE,
              data.s=ImmunoMoscatoPatients[NUMIHC,],x.label="Time from inclusion (months)",
              y.label="Probability of PFS", main.title= "Treatment with ICB \n relative to PD-L1 staining", 
              .col=c(1,2,3), o.text = "",
              .lty=c(1,1,1), show.n.risk=3, n.risk.step=1.5, n.risk.cex=0.8,
              verbose=T, frame=F)
legend("topright", lty=1, col = c(1,2,3,0,0), title = "IHC PDL1",
       legend = c("0%", "1-50%",">50%", paste("HR =", round(VALUESPFS$coefficients[2],3), "CI95[", round(VALUESPFS$conf.int[3],3), "-", round(VALUESPFS$conf.int[4],3),"]" ), 
                  paste("Logrank p =", round(VALUESPFS$sctest[3],3) )))

VALUESPFS<-summary(coxph(ImmunoMoscatoPatients$SurvOS[NUMIHC]~ImmunoMoscatoPatients$IHCPDL1[NUMIHC]))
km.coxph.plot(formula.s=Surv(ImmunoMoscatoPatients$OS[NUMIHC]/30,ImmunoMoscatoPatients$death[NUMIHC])~ImmunoMoscatoPatients$IHCPDL1[NUMIHC], 
              mark.time=TRUE,
              data.s=ImmunoMoscatoPatients[NUMIHC,],x.label="Time from inclusion (months)",
              y.label="Probability of OS", main.title= "Treatment with ICB \n relative to PD-L1 staining", 
              .col=c(1,2,3), o.text = "",
              .lty=c(1,1,1), show.n.risk=3, n.risk.step=1.5, n.risk.cex=0.8,
              verbose=T, frame=F)
legend("topright", lty=1, col = c(1,2,3,0,0), title = "IHC PDL1",
       legend = c("0%", "1-50%",">50%", paste("HR =", round(VALUESPFS$coefficients[2],3), "CI95[", round(VALUESPFS$conf.int[3],3), "-", round(VALUESPFS$conf.int[4],3),"]" ), 
                  paste("Logrank p =", round(VALUESPFS$sctest[3],3) )))


# other tests
boxplot(ImmunoMoscatoPatients$IHCPDL1[NUMIHC]~ImmunoMoscatoPatients$ORR[NUMIHC])
kruskal.test(ImmunoMoscatoPatients$ORR[NUMIHC],ImmunoMoscatoPatients$IHCPDL1[NUMIHC])
VALUEORR<-summary(glm(ImmunoMoscatoPatients$ORR[NUMIHC]~ImmunoMoscatoPatients$IHCPDL1[NUMIHC],family = quasibinomial(link = "logit")))
VALUEORR
