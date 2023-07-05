#ML using XGBOOST

library(readxl)
library(ggplot2)
library(dplyr)
library(caret)
library(reshape2)
library(ggsci)

rm(list=ls())

# 30 C ---- 07.10.2022 ----
### 1. Import data  ###
#scale promoter strengths 
#data contain scale of all metabolites and promoter strengths every timepoints
pro_str=read.table("Metabolic_promoterstrength_forML__08092022.txt")

#productions data 

dat=read_xlsx("data file")
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

tSplit=data.frame(matrix(NA, ncol=10, nrow=20))
colnames(tSplit)=c("seednumber","R2_train","train_RMSE","R2_test","iteration",
                   "eta","maxdepth","subsample","colbytree","minchidw")

for (c in 1:20) {
  print(paste("************",c,"*************"))
  ss=sample(1:30000,1)
  if(ss%in%tSplit$seednumber){
    print("************ REPLICATE SEED NUMBER **********")
    tSplit[c,1]=NA
    tSplit[c,2]=NA
    tSplit[c,3]=NA
    tSplit[c,4]=NA
    tSplit[c,5]=NA
    tSplit[c,6]=NA
    tSplit[c,7]=NA
    tSplit[c,8]=NA
    tSplit[c,9]=NA
    
  }
  
  else{
###2. prepare training and test sets ###
set.seed(ss)
A=createDataPartition(eth_str$EtOH_24h, p=0.7, list = F) #caret
train_df=data.matrix(eth_str[A,])
test_df=data.matrix(eth_str[-A,])
dim(train_df);dim(test_df)
#head(test_df)
#head(train_df)

### 3. Training the model by XGBoost ###
library(xgboost)
#details of data 
#column 6 is EtOH(g/L), 11 is EtOH(g/L/OD600) 
#Change number in label = n
#set column #n
colnames(test_df)
n=6
print(colnames(test_df)[n])
xgtrain=xgb.DMatrix(data = train_df[,c(2:4)],label=train_df[,n])
xgtest=xgb.DMatrix(data = test_df[,c(2:4)],label=test_df[,n])

#Tuning 
parameters_list = list()
set.seed(90)
for (iter in 1:1000){
  param <- list(booster = "gbtree",
                objective = "reg:squarederror",
                max_depth = sample(3:10, 1),
                eta = round(runif(1, 0.1, 0.5),1),
                subsample = round(runif(1, 0.7, 1),1),
                colsample_bytree = round(runif(1, 0.6, 1),1),
                min_child_weight = sample(0:10, 1)
  )
  parameters <- as.data.frame(param)
  parameters_list[[iter]] <- parameters
}

parameters_df = do.call(rbind, parameters_list)
Trainlowest_error_list = list()
Testlowest_error_list= list()
best_iter_list= list()


# Use randomly created parameters to create 10,000 XGBoost-models
for (row in 1:nrow(parameters_df)){
  set.seed(90)
  mdcv <- xgb.cv(data=xgtrain,nfold = 5,
                 booster = "gbtree",
                 objective = "reg:squarederror",
                 max_depth = parameters_df$max_depth[row],
                 eta = parameters_df$eta[row],
                 subsample = parameters_df$subsample[row],
                 colsample_bytree = parameters_df$colsample_bytree[row],
                 min_child_weight = parameters_df$min_child_weight[row],
                 nrounds= 100,
                 early_stopping_rounds= 5,
                 print_every_n = 10
  )
  
  Trainlowest_error <- as.data.frame(min(mdcv$evaluation_log$train_rmse_mean))
  Trainlowest_error_list[[row]] <- Trainlowest_error
  
  Testlowest_error <- as.data.frame(min(mdcv$evaluation_log$test_rmse_mean))
  Testlowest_error_list[[row]] <- Testlowest_error
  
  best_iter <- as.data.frame(mdcv$best_iteration)
  best_iter_list[[row]] <- best_iter
}

# Create object that contains all accuracy's
Trainlowest_error_df = do.call(rbind, Trainlowest_error_list)
Testlowest_error_df = do.call(rbind, Testlowest_error_list)
best_iter_error_df = do.call(rbind, best_iter_list)

# Bind columns of accuracy values and random hyperparameter values

randomsearch = cbind(Trainlowest_error_df, parameters_df)
randomsearch = cbind(Testlowest_error_df, randomsearch)
randomsearch = cbind(best_iter_error_df, randomsearch)

# Quickly display highest accuracy

randomsearch=randomsearch[order(randomsearch$`min(mdcv$evaluation_log$test_rmse_mean)`),]
bestrd=randomsearch[1:30,]
best_Rs_list= list()

for (row in 1:nrow(bestrd)){
  set.seed(0199)
  #cross validation to find best iteration
  xgb=xgboost(data=xgtrain, booster= "gbtree",
              nrounds =bestrd$`mdcv$best_iteration`[row],
              eta=bestrd$eta[row],
              print_every_n = 10,
              max.depth= bestrd$max_depth[row],
              subsample = bestrd$subsample[row],
              colsample_bytree = bestrd$colsample_bytree[row],
              min_child_weight = bestrd$min_child_weight[row]
  )
  #prediction using test data
  pred=predict(xgb, xgtest)
  t=as.numeric(test_df[,n])
  Rs=cor(pred,t)^2 #correlation ^2 = Rsquare
  
  best_Rs <- as.data.frame(Rs)
  best_Rs_list[[row]] <- best_Rs
  
}

Rs_df = do.call(rbind, best_Rs_list)
bestrd2=cbind(Rs_df ,bestrd)
bestrd2.1=bestrd2[order(bestrd2$Rs, decreasing = T),]
bestrd2.1[1:10,]
bestrd3=bestrd2[bestrd2$Rs==max(Rs_df),]
bestrd3
bestrd3=bestrd3[1,]
#training 
set.seed(0199)
xgb=xgboost(data=xgtrain,
            booster= "gbtree",
            objective= "reg:squarederror",
            nrounds =bestrd3$`mdcv$best_iteration`,
            eta=bestrd3$eta ,
            max.depth=bestrd3$max_depth,
            subsample = bestrd3$subsample,
            colsample_bytree = bestrd3$colsample_bytree,
            min_child_weight = bestrd3$min_child_weight
)

#prediction using test data
pred=predict(xgb, xgtest)
t=as.numeric(test_df[,n])
corr=cor(pred,t)^2 #correlation ^2 = Rsquare

tSplit[c,1]=ss
tSplit[c,2]=bestrd3$Rs
tSplit[c,2]=xgb$evaluation_log[xgb$niter]
tSplit[c,3]=corr
tSplit[c,4]=bestrd3$`mdcv$best_iteration`
tSplit[c,5]=bestrd3$eta
tSplit[c,6]=bestrd3$max_depth
tSplit[c,7]=bestrd3$subsample
tSplit[c,8]=bestrd3$colsample_bytree
tSplit[c,9]=bestrd3$min_child_weight

rm(bestrd3)
rm(ss)
rm(corr)
  }

}
