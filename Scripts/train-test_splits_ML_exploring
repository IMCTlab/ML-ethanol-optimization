
library(readxl)
library(ggplot2)
library(dplyr)
library(caret)
library(reshape2)
library(ggsci)

rm(list=ls())

# 30 C 
### 1. Import data  ###
#scale promoter strengths 
#data contain scale of all metabolites and promoter strengths every timepoints
pro_str=read.table("D:/BIOTEC NSTDA/ML_ethanol/Results/30C/table_scalenocenter/Metabolic_promoterstrength_forML__08092022.txt")

#productions data 
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
## 1. prepare data /spliting ##
###2. prepare training and test sets ###
# createDataPartition() for sampling in caret
tSplit=data.frame(matrix(NA, ncol=8, nrow=20))
colnames(tSplit)=c("seednumber","R2","degree","scale","C",
                  "train_RMSE","R2_train","MAE_train")

ss=sample(1:30000,20)

for(nx in 1:20){
  
  rm(r2)
  rm(svm1)
  rm(svm2)
  
  print(paste("******************",nx,"******************"))
  print(ss[nx])
  
set.seed(ss[nx])
A=createDataPartition(eth_str$EtOH_24h, p=0.7, list = F) #caret
train_df=data.matrix(eth_str[A,])
test_df=eth_str[-A,]
dim(train_df);dim(test_df)
colnames(test_df)
n=6
print(colnames(test_df)[n])
train_df2 = train_df[,c(2,3,4,6)]

## 2. Training ##
tr <- trainControl(method = "cv", number = 5)

set.seed(66)
#### SVM

svm1=train(EtOH_24h~.,train_df2, method="svmPoly", trControl= tr)
svm1$bestTune

svm2=svm1$results[svm1$results$degree%in%svm1$bestTune$degree&
               svm1$results$scale%in%svm1$bestTune$scale&
               svm1$results$C%in%svm1$bestTune$C,]

predict_svm1=as.numeric(svm1 %>% predict(test_df[,c(2,3,4)]))
preds_svm=bind_cols(Predicted=predict_svm1,Actual =test_df$EtOH_24h)
#RMSE(pred=preds_svm$Predicted,obs=preds_svm$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_svm1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_svm1)^2
mae=MAE(test_df$EtOH_24h,predict_svm1)


tSplit[nx,1]=ss[nx]
tSplit[nx,2]=r2
tSplit[nx,3]=svm2$degree
tSplit[nx,4]=svm2$scale
tSplit[nx,5]=svm2$C
tSplit[nx,6]=svm2$RMSE
tSplit[nx,7]=svm2$Rsquared
tSplit[nx,8]=svm2$MAE

}

setwd("D:/BIOTEC NSTDA/ML_ethanol/ML_ethanol_allnew/Table")
write.csv(tSplit, "Allexplore_Mutiple_train-test_split_SVM_310523.csv")


################### ML algorithms  ##############################
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

######## GLM
glm1=train(EtOH_24h~.,train_df2, method="glmnet", trControl= tr)
#glm1$results
glm2=glm1$results[glm1$results$alpha%in%glm1$bestTune$alpha&glm1$results$lambda%in%glm1$bestTune$lambda,]

predict_glm1=as.numeric(glm1 %>% predict(test_df[,c(2,3,4)]))
preds_glm=bind_cols(Predicted=predict_glm1,Actual =test_df$EtOH_24h)
RMSE(pred=preds_glm$Predicted,obs=preds_glm$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_glm1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_glm1)^2
mae=MAE(test_df$EtOH_24h,predict_glm1)

tSplit[nx,1]=ss[nx]
tSplit[nx,2]=r2
tSplit[nx,3]= glm1$bestTune$alpha
tSplit[nx,4]= glm1$bestTune$lambda
tSplit[nx,5]=glm2$RMSE
tSplit[nx,6]=glm2$Rsquared
tSplit[nx,7]=glm2$MAE

######## DT
dt1=train(EtOH_24h~.,train_df2, method="rpart", trControl= tr)
dt1$results
dt2=dt1$results[dt1$results$cp%in%dt1$bestTune$cp,]
predict_dt1=as.numeric(dt1 %>% predict(test_df[,c(2,3,4)]))
preds_dt=bind_cols(Predicted=predict_dt1,Actual =test_df$EtOH_24h)
RMSE(pred=preds_dt$Predicted,obs=preds_dt$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_dt1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_dt1)^2
mae=MAE(test_df$EtOH_24h,predict_dt1)

tSplit[nx,1]=ss[nx]
tSplit[nx,2]=r2
tSplit[nx,3]=dt1$bestTune$cp
tSplit[nx,4]=dt2$RMSE
tSplit[nx,5]=dt2$Rsquared
tSplit[nx,6]=dt2$MAE

#### Random forest
rf1=train(EtOH_24h~.,train_df2, method="rf", trControl= tr)
#rf1$bestTune
rf2=rf1$results[rf1$results$mtry%in%rf1$bestTune$mtry,]  
predict_rf1=as.numeric(rf1 %>% predict(test_df[,c(2,3,4)]))
preds_rf=bind_cols(Predicted=predict_rf1,Actual =test_df$EtOH_24h)
RMSE(pred=preds_rf$Predicted,obs=preds_rf$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_rf1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_rf1)^2
mae=MAE(test_df$EtOH_24h,predict_rf1)


tSplit[nx,1]=ss[nx]
tSplit[nx,2]=r2
tSplit[nx,3]=rf1$bestTune$mtry
tSplit[nx,4]=rf2$RMSE
tSplit[nx,5]=rf2$Rsquared
tSplit[nx,6]=rf2$MAE


#### gradient boosted tree
gb1=train(EtOH_24h~.,train_df2, method="xgbTree", trControl= tr)
#gb1$bestTune
gb2=gb1$results[rownames(gb1$results)%in%rownames(gb1$bestTune),]

predict_gb1=as.numeric(gb1 %>% predict(test_df[,c(2,3,4)]))
preds_gb=bind_cols(Predicted=predict_gb1,Actual =test_df$EtOH_24h)
#RMSE(pred=preds_gb$Predicted,obs=preds_gb$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_gb1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_gb1)^2
mae=MAE(test_df$EtOH_24h,predict_gb1)

tSplit[nx,1]=ss[nx]
tSplit[nx,2]=r2
tSplit[nx,3]=gb2$eta
tSplit[nx,4]=gb2$max_depth
tSplit[nx,5]=gb2$gamma
tSplit[nx,6]=gb2$colsample_bytree
tSplit[nx,7]=gb2$min_child_weight
tSplit[nx,8]=gb2$subsample
tSplit[nx,9]=gb2$nrounds
tSplit[nx,10]=gb2$RMSE
tSplit[nx,11]=gb2$Rsquared
tSplit[nx,12]=gb2$MAE

#### SVM
svm1=train(EtOH_24h~.,train_df2, method="svmPoly", trControl= tr)
svm1$bestTune

svm2=svm1$results[svm1$results$degree%in%svm1$bestTune$degree&
                    svm1$results$scale%in%svm1$bestTune$scale&
                    svm1$results$C%in%svm1$bestTune$C,]

predict_svm1=as.numeric(svm1 %>% predict(test_df[,c(2,3,4)]))
preds_svm=bind_cols(Predicted=predict_svm1,Actual =test_df$EtOH_24h)
#RMSE(pred=preds_svm$Predicted,obs=preds_svm$Actual)
#RMSE
rmse=sqrt(mean((test_df$EtOH_24h-predict_svm1)^2))
#R2
r2=cor(test_df$EtOH_24h,predict_svm1)^2
mae=MAE(test_df$EtOH_24h,predict_svm1)


tSplit[nx,1]=ss[nx]
tSplit[nx,2]=r2
tSplit[nx,3]=svm2$degree
tSplit[nx,4]=svm2$scale
tSplit[nx,5]=svm2$C
tSplit[nx,6]=svm2$RMSE
tSplit[nx,7]=svm2$Rsquared
tSplit[nx,8]=svm2$MAE
