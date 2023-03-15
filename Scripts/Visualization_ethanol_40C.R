
library(ggplot2)
library(reshape2)
library(readxl)
library(ggpubr)
library(dplyr)
##################

#Visualization 
# 40 C set1
setwd("D:/BIOTEC NSTDA/ML_ethanol/Results/4042C/raw")

C40=read_xlsx("Summary 40 42 oC_innova.xlsx")
C40_1=C40[-1,-c(3,5,7,9,11,13)] #cut 18h columns
colnames(C40_1)[1]="Group"
#create new dataframe with numeric 
C40_2=cbind("Name"=C40_1$Group,as.data.frame(sapply(C40_1[,2:7], as.numeric)))
C40_2$Name=substr(C40_2$Name,1,3) 
C40_2$`Glucose consumption (%)`=(C40_2$`Glucose consumption (g/L)`/100)*100
C40_2=C40_2[,c(1,2,3,8,5,6,7)]
#specific product 
C40_2$`EtOH(g/L/OD600)`=C40_2$`EtOH (g/L)`/C40_2$A600

A=C40_2%>%group_by(Name)%>%summarise_if(is.numeric,mean)
A[A$Name%in%"WT-","Name"]="xxx" #change "WT-" to "xxx" follow promoter
colnames(A)=c("Name","OD600","EtOH(g/L)","%Glucose consumption",
              "Glycerol","Pyruvate","Acetate","EtOH(g/L/OD600)")

Asd=C40_2%>%group_by(Name)%>%summarise_if(is.numeric,sd)
Asd[Asd$Name%in%"WT-","Name"]="xxx" #change "WT-" to "xxx" follow promoter
colnames(Asd)=c("Name","STD_OD600","STD_EtOH(g/L)","STD_%Glucose consumption",
                "STD_Glycerol","STD_Pyruvate","STD_Acetate","STD_EtOH(g/L/OD600)")


### 40C set2
setwd("D:/BIOTEC NSTDA/ML_ethanol/Results/40C_21sep22")

C40=read_xlsx("Summary 40 oC 24h 0.25OD starter Innova EQS 21-09-2022.xlsx")
C40_2=C40[-c(1:3),]
C40_2=C40_2[,c(1,2,3,5,6,7,8)]
colnames(C40_2)[1]="Group"
#create new dataframe with numeric 
C40_2=cbind("Name"=C40_2$Group,as.data.frame(sapply(C40_2[,2:7], as.numeric)))
C40_2$Name=substr(C40_2$Name,1,3) 
#specific product 
C40_2$`EtOH (g/L/OD600)`=C40_2[,3]/C40_2[,2]

B=C40_2%>%group_by(Name)%>%summarise_if(is.numeric,mean)
B[B$Name%in%"WT-","Name"]="xxx" #change "WT-" to "xxx" follow promoter
colnames(B)=c("Name","OD600","EtOH(g/L)","%Glucose consumption",
              "Glycerol","Pyruvate","Acetate","EtOH(g/L/OD600)")

Bsd=C40_2%>%group_by(Name)%>%summarise_if(is.numeric,sd)
Bsd[Bsd$Name%in%"WT-","Name"]="xxx" #change "WT-" to "xxx" follow promoter
colnames(Bsd)=c("Name","STD_OD600","STD_EtOH(g/L)","STD_%Glucose consumption",
                "STD_Glycerol","STD_Pyruvate","STD_Acetate","STD_EtOH(g/L/OD600)")


C=rbind(A,B)
C=C[-10,]  #remove xxx(WT) from setA
Csd=rbind(Asd,Bsd)
Csd=Csd[-10,]  #remove xxx(WT) from setA

All=cbind(C,Csd)
All1=All[,c(1,2,10,3,11,4,12,5,13,6,14,7,15,8,16)]
setwd("D:/BIOTEC NSTDA/ML_ethanol/Manuscript2023/Figure/Table")
writexl::write_xlsx(All1, path =  "Rawdata_40CEthanolandMetabolites.xlsx")

####  plot2 #EtOH (g/L/OD600) vs other 
p1=ggplot(C, aes(`EtOH(g/L/OD600)`, OD600, color=Name, label=Name))+
  geom_point(size=2.5)+
  ggtitle("40C EtOH production vs OD600")+
  xlim(0,round(max(C$`EtOH(g/L/OD600)`)+0.5))+ylim(0,round(max(C$OD600))+0.5)+
  geom_text(size=2,hjust =0.5, vjust=1.5, color="black")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

p2=ggplot(C, aes(`EtOH(g/L/OD600)`,`%Glucose consumption`, color=Name, label=Name))+
  geom_point(size=2.5, )+
  ggtitle("40C EtOH production vs %Glucose consumption")+
  xlim(0,round(max(C$`EtOH(g/L/OD600)`)+0.5))+ylim(0,round(max(C$`%Glucose consumption`))+0.5)+
  geom_text(size=2,hjust =0.5, vjust=1.5, color="black")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

p3=ggplot(C, aes(`EtOH(g/L/OD600)`, `EtOH(g/L)`, color=Name, label=Name))+
  geom_point(size=2.5)+
  ggtitle("40C EtOH (g/L/OD600) vs EtOH (g/L)")+
  xlim(0,round(max(C$`EtOH(g/L/OD600)`)+0.5))+ylim(0,round(max(C$`EtOH(g/L)`))+0.5)+
  geom_text(size=2,hjust =0.5, vjust=1.5, color="black")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

p4=ggplot(C, aes(`EtOH(g/L/OD600)`, Glycerol, color=Name, label=Name))+
  geom_point(size=2.5)+
  ggtitle("40C EtOH production vs Glycerol")+
  xlim(0,round(max(C$`EtOH(g/L/OD600)`)+0.5))+ylim(0,round(max(C$Glycerol))+0.5)+
  geom_text(size=2,hjust =0.5, vjust=1.5, color="black")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

p5=ggplot(C, aes(`EtOH(g/L/OD600)`, Pyruvate, color=Name, label=Name))+
  geom_point(size=2.5)+
  ggtitle("40C EtOH production vs Pyruvate")+
  xlim(0,round(max(C$`EtOH(g/L/OD600)`)+0.5))+ylim(0,round(max(C$Pyruvate))+0.5)+
  geom_text(size=2,hjust =0.5, vjust=1.5, color="black")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

p6=ggplot(C, aes(`EtOH(g/L/OD600)`, Acetate, color=Name, label=Name))+
  geom_point(size=2.5)+
  ggtitle("40C EtOH production vs Acetate")+
  xlim(0,round(max(C$`EtOH(g/L/OD600)`)+0.5))+ylim(0,round(max(C$Acetate))+0.5)+
  geom_text(size=2,hjust =0.5, vjust=1.5, color="black")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")


p_eth40C=ggarrange(p3,p1,p2,p4,p5,p6,ncol=3, nrow = 2)
P2=annotate_figure(p_eth40C, top = text_grob("production at 40C, 24h ", 
                                             color = "black", face = "bold", size = 14))
