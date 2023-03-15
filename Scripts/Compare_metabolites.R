
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
colnames(Asd)=c("Name","OD600","EtOH(g/L)","%Glucose consumption",
                "Glycerol","Pyruvate","Acetate","EtOH(g/L/OD600)")


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
colnames(Bsd)=c("Name","OD600","EtOH(g/L)","%Glucose consumption",
                "Glycerol","Pyruvate","Acetate","EtOH(g/L/OD600)")


C=rbind(A,B)
C=C[-10,]  #remove xxx(WT) from setA


## Add 30 C data ##

#Visualization 
#rep1 30 C
setwd("D:/BIOTEC NSTDA/ML_ethanol/Results/30C/raw/")
dat=read_xlsx("Metabolite and growth data_Day 1 and 2.xlsx", sheet = 1)
#define column name
colnames(dat)=c("Name","A600_24h","A600_24h_std","A600_48h","A600_48h_std",
                "EtOH_24h","EtOH_24h_std","EtOH_48h","EtOH_48h_std",
                "Glu_24h","Glu_24h_std","Glu_48h","Glu_48h_std",
                "Gly_24h","Gly_24h_std","Gly_48h","Gly_48h_std",
                "Acetate_24h","Acetate_24h_std","Acetate_48h","Acetate_48h_std",
                "Pyruvate_24h","Pyruvate_24h_std","Pyruvate_48h","Pyruvate_48h_std")
dat=dat[-c(1,2),] #remove first two rows
dat1=as.data.frame(sapply(dat[,2:25], as.numeric)) #change to numeric
dat1$Name=dat$Name #copy rownames to new "dat1"
rownames(dat1)=dat1$Name

mdatsd=melt(dat1[,c(colnames(dat1)[grepl("std",colnames(dat1))],"Name")])
dat1=dat1[,colnames(dat1)[!grepl("std",colnames(dat1))]]
dat1_24=dat1[,grepl("24h",colnames(dat1))]
dat1_24$Name=row.names(dat1_24)
dat1_24$`EtOH_24h(g/L/OD600)`= dat1_24$EtOH_24h/dat1_24$A600_24h

#####
#Ethanol
et=left_join(C[,c(1,3)],dat1_24[,c(2,7)])
colnames(et)=c("Name","40C","30C")
et=et[order(et$`40C`, decreasing = T),]
X=factor(et$Name, levels=et$Name)
et$Name=X
met=melt(et)
p1=ggplot(met, aes(Name, value, color=variable,group=variable))+
  geom_line(size=1.5)+geom_point(size=2)+
  ylab("Ethanol (g/L)")+labs(color="conditions")+
  theme_bw()+ggtitle("Ethanol")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())

#Specific product
et=left_join(C[,c(1,8)],dat1_24[,c(8,7)])
colnames(et)=c("Name","40C","30C")
et=et[order(et$`40C`, decreasing = T),]
et$Name=factor(et$Name, levels=levels(X))
met=melt(et)
p2=ggplot(met, aes(Name, value, color=variable,group=variable))+
  geom_line(size=1.5)+geom_point(size=2)+
  ylab(expression(paste("Ethanol (g/L/",OD[600],")")))+labs(color="conditions")+
  theme_bw()+ggtitle("Specfic product of ethanol")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())

#Pyruvate
et=left_join(C[,c(1,6)],dat1_24[,c(7,6)])
colnames(et)=c("Name","40C","30C")
et=et[order(et$`40C`, decreasing = T),]
et$Name=factor(et$Name, levels=levels(X))
met=melt(et)
p3=ggplot(met, aes(Name, value, color=variable,group=variable))+
  geom_line(size=1.5)+geom_point(size=2)+
  ylab("Pyruvate (g/L)")+labs(color="conditions")+
  theme_bw()+ggtitle("Pyruvate")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size=15), 
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())


#Acetate
et=left_join(C[,c(1,7)],dat1_24[,c(7,5)])
colnames(et)=c("Name","40C","30C")
et=et[order(et$`40C`, decreasing = T),]
et$Name=factor(et$Name, levels=levels(X))
met=melt(et)
p4=ggplot(met, aes(Name, value, color=variable,group=variable))+
  geom_line(size=1.5)+geom_point(size=2)+
  ylab("Acetate (g/L)")+labs(color="conditions")+
  theme_bw()+ggtitle("Acetate")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())


#Glycerol
et=left_join(C[,c(1,5)],dat1_24[,c(7,4)])
colnames(et)=c("Name","40C","30C")
et=et[order(et$`40C`, decreasing = T),]
et$Name=factor(et$Name, levels=levels(X))
met=melt(et)
p5=ggplot(met, aes(Name, value, color=variable,group=variable))+
  geom_line(size=1.5)+geom_point(size=2)+
  ylab("Glycerol (g/L)")+labs(color="conditions")+
  theme_bw()+ggtitle("Glycerol ")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())

#Glucose consumption
et=left_join(C[,c(1,4)],dat1_24[,c(7,3)])
colnames(et)=c("Name","40C","30C")
et=et[order(et$`40C`, decreasing = T),]
et$Name=factor(et$Name, levels=levels(X))
met=melt(et)
p6=ggplot(met, aes(Name, value, color=variable,group=variable))+
  geom_line(size=1.5)+geom_point(size=2)+
  ylab("% Glucose consumption")+labs(color="conditions")+
  theme_bw()+ggtitle("% Glucose consumption")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())

#OD600
et=left_join(C[,c(1,2)],dat1_24[,c(7,1)])
colnames(et)=c("Name","40C","30C")
et=et[order(et$`40C`, decreasing = T),]
et$Name=factor(et$Name, levels=levels(X))
met=melt(et)
p7=ggplot(met, aes(Name, value, color=variable,group=variable))+
  geom_line(size=1.5)+geom_point(size=2)+
  ylab(bquote(OD[600]))+labs(color="conditions")+
  theme_bw()+ggtitle(bquote(OD[600]))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())

pp=ggarrange(p1,p2,p7,p6,p3,p4,p5,ncol=3, nrow = 3)
pp
