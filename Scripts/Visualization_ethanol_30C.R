
library(ggplot2)
library(reshape2)
library(readxl)
library(ggpubr)

##################

#Visualization 
#rep1 30 C

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



####  plot2 #EtOH vs other 
c=cor.test(dat1_24$EtOH_24h,dat1_24$A600_24h)
p1=ggplot(dat1_24, aes(EtOH_24h, A600_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  annotate(geom="text", size=6,x=16,y=32,label=paste("Pearson's correlation=",round(c$estimate,2)))+
  #ggtitle("30C EtOH production vs OD600")+
  xlab("Ethanol production (g/L)")+ylab(bquote(OD[600]))+
  xlim(0,round(max(dat1_24$EtOH_24h)+0.5))+ylim(0,round(max(dat1_24$A600_24h)+0.5))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        text = element_text(size=20))

c=cor.test(dat1_24$EtOH_24h,dat1_24$Glu_24h)
p2=ggplot(dat1_24, aes(EtOH_24h, Glu_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7 )+
  annotate(geom="text",size=6,x=16,y=96,label=paste("Pearson's correlation=",round(c$estimate,2)))+
  #ggtitle("30C EtOH production vs %Glu consumption")+
  xlab("Ethanol production (g/L)")+ylab("% Glucose comsumption")+
  xlim(0,round(max(dat1_24$EtOH_24h)+0.5))+ylim(0,round(max(dat1_24$Glu_24h)+0.5))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        text = element_text(size=20))

c=cor.test(dat1_24$EtOH_24h,dat1_24$`EtOH_24h(g/L/OD600)`)
p3=ggplot(dat1_24, aes(EtOH_24h, `EtOH_24h(g/L/OD600)`, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  annotate(geom="text", size=6,x=16,y=2.8,label=paste("Pearson's correlation=",round(c$estimate,2)))+
  #ggtitle("30C EtOH (g/L/OD600) vs EthOH (g/L)")+
  xlab("Ethanol production (g/L)")+
  ylab(expression(paste("Specific product of ethanol (g/L/",OD[600],")")))+
  xlim(0,round(max(dat1_24$EtOH_24h)+0.5))+ylim(0,round(max(dat1_24$`EtOH_24h(g/L/OD600)`)+0.5))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        text = element_text(size=20))

c=cor.test(dat1_24$EtOH_24h,dat1_24$Gly_24h)
p4=ggplot(dat1_24, aes(EtOH_24h, Gly_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  annotate(geom="text", size=6,x=16,y=4.6,label=paste("Pearson's correlation=",round(c$estimate,2)))+
  #ggtitle("30C EtOH production vs Glycerol")+
  xlab("Ethanol production (g/L)")+ylab("Glycerol (g/L)")+
  xlim(0,round(max(dat1_24$EtOH_24h)+0.5))+ylim(0,round(max(dat1_24$Gly_24h)+0.5))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        text = element_text(size=20))

c=cor.test(dat1_24$EtOH_24h,dat1_24$Pyruvate_24h)
p5=ggplot(dat1_24, aes(EtOH_24h, Pyruvate_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  annotate(geom="text", size=6,x=16,y=1.85,label=paste("Pearson's correlation=",round(c$estimate,2)))+
  #ggtitle("30C EtOH production vs Pyruvate")+
  xlab("Ethanol production (g/L)")+ylab("Pyruvate (g/L)")+
  xlim(0,round(max(dat1_24$EtOH_24h)+0.5))+ylim(0,round(max(dat1_24$Pyruvate)+0.5))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        text = element_text(size=20))

c=cor.test(dat1_24$EtOH_24h,dat1_24$Acetate_24h)
p6=ggplot(dat1_24, aes(EtOH_24h, Acetate_24h, color=Name, label=Name))+
  geom_point(size=2.5, alpha=0.7)+
  annotate(geom="text", size=6,x=16,y=2.8,label=paste("Pearson's correlation=",round(c$estimate,2)))+
  #ggtitle("30C EtOH production vs Acetate")+
  xlab("Ethanol production (g/L)")+ylab("Acetate (g/L)")+
  xlim(0,round(max(dat1_24$EtOH_24h)+0.5))+ylim(0,round(max(dat1_24$Acetate_24h)+0.5))+
  #geom_text(size=2,hjust =0.5, vjust=1.5, color="grey")
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        text = element_text(size=20))


p_eth40C=ggarrange(p3,p1,p2,p4,p5,p6,ncol=3, nrow = 2)
P2=annotate_figure(p_eth40C)
                          
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
        axis.ticks.x= element_blank(),
        text = element_text(size=20))+
  geom_hline(yintercept = dat1_24_2[dat1_24_2$Name%in%"xxx","Pyruvate_24h"],
              linetype="dashed")

w2=ggplot(dat1_24_2, aes(x=Name,y=Acetate_24h, group=1))+
  geom_line(size=1)+ #, color="#F8766D"
  xlab("strains")+ylab("Acetate (g/L)")+
  scale_y_continuous(breaks = seq(0,3,0.5), limits = c(0,3))+
  ggtitle("Acetate")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        text = element_text(size=20))+
  geom_hline(yintercept =dat1_24_2[dat1_24_2$Name%in%"xxx","Acetate_24h"],
              linetype="dashed")

w3=ggplot(dat1_24_2, aes(x=Name,y=Gly_24h, group=1))+
  geom_line(size=1)+ #, color="#7CAE00"
  xlab("strains")+ylab("Glycerol (g/L)")+
  scale_y_continuous(breaks = seq(0,3,0.5), limits = c(0,3))+
  ggtitle("Glycerol")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        text = element_text(size=20))+
  geom_hline(yintercept = dat1_24_2[dat1_24_2$Name%in%"xxx","Gly_24h"],
              linetype="dashed")

w4=ggplot(dat1_24_2, aes(x=Name,y=`EtOH_24h(g/L/OD600)`, group=1))+
  geom_line(size=1)+ #, color="#C77CFF"
  xlab("strains")+ylab(expression(paste("ethanol (g/L/",OD[600],")")))+
  scale_y_continuous(breaks = seq(0,3,0.5), limits = c(0,3))+
  ggtitle("Specific product of ethanol")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        text = element_text(size=20))+
  geom_hline(yintercept =  dat1_24_2[dat1_24_2$Name%in%"xxx","EtOH_24h(g/L/OD600)"], 
              linetype="dashed")


w5=ggplot(dat1_24_2, aes(x=Name,y=EtOH_24h, group=1))+
  geom_line(size=1)+ #, color="#E68613"
  xlab("strains")+ylab("Ethanol production (g/L)")+
  scale_y_continuous(breaks = seq(0,70,10), limits = c(0,70))+
  ggtitle("Ethanol")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        text = element_text(size=20))+
  geom_hline(yintercept =  dat1_24_2[dat1_24_2$Name%in%"xxx","EtOH_24h"], 
             linetype="dashed")

w6=ggplot(dat1_24_2, aes(x=Name,y=Glu_24h, group=1))+
  geom_line(size=1)+ #, color="#00A9FF"
  xlab("strains")+ylab("%Glucose consumption")+
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100))+
  ggtitle("%Glucose consumption")+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        text = element_text(size=20))+
  geom_hline(yintercept =  dat1_24_2[dat1_24_2$Name%in%"xxx","Glu_24h"],
              linetype="dashed")

w7=ggplot(dat1_24_2, aes(x=Name,y=A600_24h, group=1))+
  geom_line(size=1)+ #, color="#FF61CC"
  xlab("strains")+ylab(bquote(OD[600]))+
  scale_y_continuous(breaks = seq(0,70,10), limits = c(0,70))+
  ggtitle(bquote(OD[600]))+theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        text = element_text(size=20))+
  geom_hline(yintercept = dat1_24_2[dat1_24_2$Name%in%"xxx","A600_24h"], 
              linetype="dashed")

ww=ggarrange(w5,w4,w6,w7,w1,w2,w3,ncol=3, nrow = 3)
