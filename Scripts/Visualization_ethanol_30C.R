
library(ggplot2)
library(reshape2)
library(readxl)
library(ggpubr)
library(ggsci)
##################

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

X=dat1_24[dat1_24$EtOH_24h>37.83359,] #xxx
X2=X[order(X$EtOH_24h, decreasing = T),]
log2(61.96128/37.83359)
log2(52.86604/37.83359)
####  plot1 #
P1=ggpairs(dat1_24[,-7])  #or ggpairs(dat1_24, columns=1:6)

####  plot2 #EtOH vs other 
p1=ggplot(dat1_24, aes(`EtOH_24h(g/L/OD600)`, A600_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  #ggtitle("30C EtOH production vs OD600")+
  xlab("Specific product of ethanol (g/L/OD600)")+ylab("OD600")+
  xlim(0,round(max(dat1_24$`EtOH_24h(g/L/OD600)`)+0.5))+ylim(0,round(max(dat1_24$A600_24h)))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")
#p1+scale_color_jco()
p2=ggplot(dat1_24, aes(`EtOH_24h(g/L/OD600)`, Glu_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7 )+
  #ggtitle("30C EtOH production vs %Glu consumption")+
  xlab("Specific product of ethanol (g/L/OD600)")+ylab("% Glucose comsumption")+
  xlim(0,round(max(dat1_24$`EtOH_24h(g/L/OD600)`)+0.5))+ylim(0,round(max(dat1_24$Glu_24h)))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

p3=ggplot(dat1_24, aes(`EtOH_24h(g/L/OD600)`, EtOH_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  #ggtitle("30C EtOH (g/L/OD600) vs EthOH (g/L)")+
  xlab("Specific product of ethanol (g/L/OD600)")+ylab("Ethanol production (g/L)")+
  xlim(0,round(max(dat1_24$`EtOH_24h(g/L/OD600)`)+0.5))+ylim(0,round(max(dat1_24$EtOH_24h)))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

p4=ggplot(dat1_24, aes(`EtOH_24h(g/L/OD600)`, Gly_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  #ggtitle("30C EtOH production vs Glycerol")+
  xlab("Specific product of ethanol (g/L/OD600)")+ylab("Glycerol (g/L)")+
  xlim(0,round(max(dat1_24$`EtOH_24h(g/L/OD600)`)+0.5))+ylim(0,round(max(dat1_24$Gly_24h)))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

p5=ggplot(dat1_24, aes(`EtOH_24h(g/L/OD600)`, Pyruvate_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  #ggtitle("30C EtOH production vs Pyruvate")+
  xlab("Specific product of ethanol (g/L/OD600)")+ylab("Pyruvate (g/L)")+
  xlim(0,round(max(dat1_24$`EtOH_24h(g/L/OD600)`)+0.5))+ylim(0,round(max(dat1_24$Pyruvate)))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

p6=ggplot(dat1_24, aes(`EtOH_24h(g/L/OD600)`, Acetate_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  #ggtitle("30C EtOH production vs Acetate")+
  xlab("Specific product of ethanol (g/L/OD600)")+ylab("Acetate (g/L)")+
  xlim(0,round(max(dat1_24$`EtOH_24h(g/L/OD600)`)+0.5))+ylim(0,round(max(dat1_24$Acetate_24h)))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")


p_eth40C=ggarrange(p3,p1,p2,p4,p5,p6,ncol=3, nrow = 2)
P2=annotate_figure(p_eth40C, top = text_grob("production at 30 C, 24h ", 
                                          color = "black", face = "bold", size = 14))
                                          
                                          
# Separate plots
n=as.data.frame(dat1_24[order(dat1_24$EtOH_24h, decreasing = T),"Name"])
Xn=as.data.frame(dat1_24[order(dat1_24$EtOH_24h, decreasing = T),])
dat1_24$Name=factor(dat1_24$Name, levels =n[,1] )    
dat1_24_2=dat1_24
w1=ggplot(dat1_24_2, aes(x=Name,y=Pyruvate_24h, group=1))+
  geom_line(size=1)+ #color="#00B8E7"
  xlab("strains")+ylab("Pyruvate (g/L)")+
  scale_y_continuous(breaks = seq(0,3,0.5), limits = c(0,3))+
  ggtitle("Pyruvate")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())+
  geom_hline(yintercept = dat1_24_2[dat1_24_2$Name%in%"xxx","Pyruvate_24h"],
              linetype="dashed")

w2=ggplot(dat1_24_2, aes(x=Name,y=Acetate_24h, group=1))+
  geom_line(size=1)+ #, color="#F8766D"
  xlab("strains")+ylab("Acetate (g/L)")+
  scale_y_continuous(breaks = seq(0,3,0.5), limits = c(0,3))+
  ggtitle("Acetate")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())+
  geom_hline(yintercept =dat1_24_2[dat1_24_2$Name%in%"xxx","Acetate_24h"],
              linetype="dashed")

w3=ggplot(dat1_24_2, aes(x=Name,y=Gly_24h, group=1))+
  geom_line(size=1)+ #, color="#7CAE00"
  xlab("strains")+ylab("Glycerol (g/L)")+
  scale_y_continuous(breaks = seq(0,3,0.5), limits = c(0,3))+
  ggtitle("Glycerol")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())+
  geom_hline(yintercept = dat1_24_2[dat1_24_2$Name%in%"xxx","Gly_24h"],
              linetype="dashed")

w4=ggplot(dat1_24_2, aes(x=Name,y=`EtOH_24h(g/L/OD600)`, group=1))+
  geom_line(size=1)+ #, color="#C77CFF"
  xlab("strains")+ylab(expression(paste("ethanol (g/L/",OD[600],")")))+
  scale_y_continuous(breaks = seq(0,3,0.5), limits = c(0,3))+
  ggtitle("Specific product of ethanol)")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())+
  geom_hline(yintercept =  dat1_24_2[dat1_24_2$Name%in%"xxx","EtOH_24h(g/L/OD600)"], 
              linetype="dashed")


w5=ggplot(dat1_24_2, aes(x=Name,y=EtOH_24h, group=1))+
  geom_line(size=1)+ #, color="#E68613"
  xlab("strains")+ylab("Ethanol production (g/L)")+
  scale_y_continuous(breaks = seq(0,70,10), limits = c(0,70))+
  ggtitle("Ethanol")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())+
  geom_hline(yintercept =  dat1_24_2[dat1_24_2$Name%in%"xxx","EtOH_24h"], 
             linetype="dashed")

w6=ggplot(dat1_24_2, aes(x=Name,y=Glu_24h, group=1))+
  geom_line(size=1)+ #, color="#00A9FF"
  xlab("strains")+ylab("%Glucose consumption")+
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100))+
  ggtitle("%Glucose consumption")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())+
  geom_hline(yintercept =  dat1_24_2[dat1_24_2$Name%in%"xxx","Glu_24h"],
              linetype="dashed")

w7=ggplot(dat1_24_2, aes(x=Name,y=A600_24h, group=1))+
  geom_line(size=1)+ #, color="#FF61CC"
  xlab("strains")+ylab(bquote(OD[600]))+
  scale_y_continuous(breaks = seq(0,70,10), limits = c(0,70))+
  ggtitle(bquote(OD[600]))+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank())+
  geom_hline(yintercept = dat1_24_2[dat1_24_2$Name%in%"xxx","A600_24h"], 
              linetype="dashed")

ww=ggarrange(w5,w4,w6,w7,w1,w2,w3,ncol=3, nrow = 3)
ww


