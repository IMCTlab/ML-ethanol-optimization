
library(readxl)
library(ggplot2)
library(dplyr)
library(caret)
library(reshape2)

# 30 C 
### 1. Import data  ###
#scale promoter strengths 
#data contain scale of all metabolites and promoter strengths every timepoints
pro_str=read.table("promoter strength file.txt")

#productions data 
setwd("path")
dat=read_xlsx("production data.xlsx", sheet = 1)
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
#dat1 contain all data and std of metabolites

mdatsd=melt(dat1[,c(colnames(dat1)[grepl("std",colnames(dat1))],"Name")])
dat1=dat1[,colnames(dat1)[!grepl("std",colnames(dat1))]]  #select average data
dat1_24=dat1[,grepl("24h",colnames(dat1))] # select time 24h
dat1_24$Name=row.names(dat1_24)
dat1_24$`EtOH_24h(g/L/OD600)`= dat1_24$EtOH_24h/dat1_24$A600_24h
# at this point dat1_24 contains all data from 24h in numeric class
pro_str2=pro_str[,c(13,23,24,25)] #select only name and promoter at 24h
eth_str=left_join(pro_str2, dat1_24, by="Name")

### ML ####
library(caret)

###1. prepare training and test sets ###
# createDataPartition() for sampling in caret
set.seed(019)
A=createDataPartition(eth_str$EtOH_24h, p=0.7, list = F) #caret
train_df=data.matrix(eth_str[A,])
test_df=eth_str[-A,]
dim(train_df);dim(test_df)

colnames(test_df)
n=6
print(colnames(test_df)[n])
train_df2=train_df[,c(2,3,4,6)]

## 2. Training ##
tr <- trainControl(method = "cv", number = 5)

#### linear regression
set.seed(66)
lm1=train(EtOH_24h~.,train_df2, method="lm", trControl= tr)
lm1$results

predict_lm1=as.numeric(lm1 %>% predict(test_df[,c(2,3,4)]))
preds_lm=bind_cols(Predicted=predict_lm1,Actual =test_df$EtOH_24h)
RMSE(pred=preds_lm$Predicted,obs=preds_lm$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_lm1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_lm1)^2
mae=MAE(test_df$EtOH_24h,predict_lm1)

p1=preds_lm%>%ggplot(aes(x=Actual, y=Predicted))+   
  geom_point(alpha = 0.8, color = "cadetblue", size=2) +
  geom_smooth(method = "loess", formula = "y ~ x", se=F, color="black") +
  coord_fixed()+xlab("Actual ethanol (g/L)")+ylab("predicted ethanol (g/L)")+
  scale_x_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  scale_y_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  theme_bw()+
  #geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(title = "Linear regression",
  subtitle= substitute(paste(R^2, " = ",r2a,
                             ", RMSE = ",rmsea), 
                       list(r2a=round(r2,4),rmsea=round(rmse,4))))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15))

#### Glm
set.seed(66)
glm1=train(EtOH_24h~.,train_df2, method="glmnet", trControl= tr)
glm1$results

predict_glm1=as.numeric(glm1 %>% predict(test_df[,c(2,3,4)]))
preds_glm=bind_cols(Predicted=predict_glm1,Actual =test_df$EtOH_24h)
RMSE(pred=preds_glm$Predicted,obs=preds_glm$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_glm1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_glm1)^2
mae=MAE(test_df$EtOH_24h,predict_glm1)

p2=preds_glm%>%ggplot(aes(x=Actual, y=Predicted))+   
  geom_point(alpha = 0.8, color = "cadetblue", size=2) +
  geom_smooth(method = "loess", formula = "y ~ x", se=F, color="black") +
  coord_fixed()+xlab("Actual ethanol (g/L)")+ylab("predicted ethanol (g/L)")+
  scale_x_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  scale_y_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  theme_bw()+
  #geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(title = " Generalized linear regression",
       subtitle= substitute(paste(R^2, " = ",r2a,
                                  ", RMSE = ",rmsea), 
                            list(r2a=round(r2,4),rmsea=round(rmse,4))))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15))

#### Decision tree
set.seed(66)
dt1=train(EtOH_24h~.,train_df2, method="rpart", trControl= tr)
dt1$results

predict_dt1=as.numeric(dt1 %>% predict(test_df[,c(2,3,4)]))
preds_dt=bind_cols(Predicted=predict_dt1,Actual =test_df$EtOH_24h)
RMSE(pred=preds_dt$Predicted,obs=preds_dt$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_dt1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_dt1)^2
mae=MAE(test_df$EtOH_24h,predict_dt1)


p3=preds_dt%>%ggplot(aes(x=Actual, y=Predicted))+   
  geom_point(alpha = 0.8, color = "cadetblue", size=2) +
  geom_smooth(method = "loess", formula = "y ~ x", se=F, color="black") +
  coord_fixed()+xlab("Actual ethanol (g/L)")+ylab("predicted ethanol (g/L)")+
  scale_x_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  scale_y_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  theme_bw()+
  #geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(title = " Decision tree",
       subtitle= substitute(paste(R^2, " = ",r2a,
                                  ", RMSE = ",rmsea), 
                            list(r2a=round(r2,4),rmsea=round(rmse,4))))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15))

#### Random forest
set.seed(66)
rf1=train(EtOH_24h~.,train_df2, method="rf", trControl= tr)
rf1$results

predict_rf1=as.numeric(rf1 %>% predict(test_df[,c(2,3,4)]))
preds_rf=bind_cols(Predicted=predict_rf1,Actual =test_df$EtOH_24h)
RMSE(pred=preds_rf$Predicted,obs=preds_rf$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_rf1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_rf1)^2
mae=MAE(test_df$EtOH_24h,predict_rf1)


p4=preds_rf%>%ggplot(aes(x=Actual, y=Predicted))+   
  geom_point(alpha = 0.8, color = "cadetblue", size=2) +
  geom_smooth(method = "loess", formula = "y ~ x", se=F, color="black") +
  coord_fixed()+xlab("Actual ethanol (g/L)")+ylab("predicted ethanol (g/L)")+
  scale_x_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  scale_y_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  theme_bw()+
  #geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(title = " Random Forest",
       subtitle= substitute(paste(R^2, " = ",r2a,
                                  ", RMSE = ",rmsea), 
                            list(r2a=round(r2,4),rmsea=round(rmse,4))))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15))


#### gradient boosted tree
set.seed(66)
gb1=train(EtOH_24h~.,train_df2, method="xgbTree", trControl= tr)
gb1$results

predict_gb1=as.numeric(gb1 %>% predict(test_df[,c(2,3,4)]))
preds_gb=bind_cols(Predicted=predict_gb1,Actual =test_df$EtOH_24h)
RMSE(pred=preds_gb$Predicted,obs=preds_gb$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_gb1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_gb1)^2
mae=MAE(test_df$EtOH_24h,predict_gb1)


p5=preds_gb%>%ggplot(aes(x=Actual, y=Predicted))+   
  geom_point(alpha = 0.8, color = "cadetblue", size=2) +
  geom_smooth(method = "loess", formula = "y ~ x", se=F, color="black") +
  coord_fixed()+xlab("Actual ethanol (g/L)")+ylab("predicted ethanol (g/L)")+
  scale_x_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  scale_y_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  theme_bw()+
  #geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(title = "Gradient boosted tree",
  subtitle= substitute(paste(R^2, " = ",r2a,
                             ", RMSE = ",rmsea), 
                       list(r2a=round(r2,4),rmsea=round(rmse,4))))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15))



#### SVM
set.seed(66)
svm1=train(EtOH_24h~.,train_df2, method="svmPoly", trControl= tr)
svm1$results

predict_svm1=as.numeric(svm1 %>% predict(test_df[,c(2,3,4)]))
preds_svm=bind_cols(Predicted=predict_svm1,Actual =test_df$EtOH_24h)
RMSE(pred=preds_svm$Predicted,obs=preds_svm$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_svm1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_svm1)^2
mae=MAE(test_df$EtOH_24h,predict_svm1)


p6=preds_svm%>%ggplot(aes(x=Actual, y=Predicted))+   
  geom_point(alpha = 0.6, color = "cadetblue", size=2) +
  geom_smooth(method = "loess", formula = "y ~ x", se=F, color="black") +
  coord_fixed()+xlab("Actual ethanol (g/L)")+ylab("predicted ethanol (g/L)")+
  scale_x_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  scale_y_continuous(limits=c(0,65),breaks = seq(0,65,by=20))+
  theme_bw()+
  #geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(title = "Support vector machine",
       subtitle= substitute(paste(R^2, " = ",r2a,
                                  ", RMSE = ",rmsea), 
                            list(r2a=round(r2,4),rmsea=round(rmse,4))))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15))

library(ggpubr)
pp=ggarrange(p1,p2,p3,p4,p5,p6,ncol=3, nrow = 2)
pp
