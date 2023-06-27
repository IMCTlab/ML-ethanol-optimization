#ML using XGBOOST

library(readxl)
library(ggplot2)
library(dplyr)
library(caret)
library(reshape2)

rm(list=ls())
### High temperature ###

### 1. Import data  ###
#scale promoter strengths 
#40C
prost2_40=read.table("Promoterstrength_40C_vectornorm_19072022.txt")
rownames(prost2_40)=c("vector","PDC1","ADH1","TPS1","ACT1",
                      "PGK1","ENO2","TDH3","YEF3")
sprost40=as.data.frame(scale(prost2_40[-1,7], center = F, scale = T))
sprost40=cbind(rownames(prost2_40)[-1],sprost40)
colnames(sprost40)=c("promoter","promoterstrength")

#productions data 
# 40 C set1
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
C=C[-c(10,27),]  #remove xxx(WT) from setA and B


# 40 C set3 validation
#setwd("D:/BIOTEC NSTDA/ML_ethanol/Results/4042C_veridation_29dec22/")
#V40=read_xlsx("Summary Predicted 15 clones 40 oC 24h (3rd flr EQS).xlsx")

V40=read_xlsx("Summary 15 Predicted Strains 40 oC 0.25OD starter 24h Innova EQS 2nd flr (3rd repeat)06-01-2023.xlsx")

V40_1=V40[-c(1,2,3,4),-5] #cut 18h columns
colnames(V40_1)[1]="Group"
#create new dataframe with numeric 
V40_2=cbind("Name"=V40_1$Group,as.data.frame(sapply(V40_1[,2:8], as.numeric)))
V40_2$Name=substr(V40_2$Name,1,3) 
#C40_2$`Glucose consumption (%)`=(C40_2$`Glucose consumption (g/L)`/100)*100
V40_2=V40_2[,c(1,2,3,5,6,8,7,4)]
#specific product 
#C40_2$`EtOH(g/L/OD600)`=C40_2$`EtOH (g/L)`/C40_2$A600

V=V40_2%>%group_by(Name)%>%summarise_if(is.numeric,mean)
V[V$Name%in%"WT-","Name"]="xxx" #change "WT-" to "xxx" follow promoter
colnames(V)=c("Name","OD600","EtOH(g/L)","%Glucose consumption",
              "Glycerol","Pyruvate","Acetate","EtOH(g/L/OD600)")

Vsd=V40_2%>%group_by(Name)%>%summarise_if(is.numeric,sd)
Vsd[Vsd$Name%in%"WT-","Name"]="xxx" #change "WT-" to "xxx" follow promoter
colnames(Vsd)=c("Name","STD_OD600","STD_EtOH(g/L)","STD_%Glucose consumption",
                "STD_Glycerol","STD_Pyruvate","STD_Acetate","STD_EtOH(g/L/OD600)")

xxx=rbind(A[A$Name=="xxx",],B[B$Name=="xxx",],V[V$Name=="xxx",])
xxx2=xxx%>%group_by(Name)%>%summarise_if(is.numeric,mean)

C=rbind(xxx2,C,V[-which(V$Name=="xxx"),])
allV=cbind(V,Vsd)
allV2=allV[,c(1,2,10,3,11,4,12,5,13,6,14,7,15,8,16)]
#setwd("D:/BIOTEC NSTDA/ML_ethanol/ML_ethanol_allnew/Table")
#writexl::write_xlsx(allV2, path = "Validation40C_averageandSTD020323.xlsx")
#add promoter strength
#import label from P'Wut
setwd("D:/BIOTEC NSTDA/ML_ethanol/Results/30C")
lb=read.csv("input_AI_Guided.csv")
colnames(lb)[1]="Name"
lb1=lb[,c(1,3,4,5)]
col=list(c("PDC1","ADH1","TPS1"))
lb1[,col[[1]]]=sprost40$promoterstrength[match(unlist(lb1[,col[[1]]]),sprost40$promoter)]

eth_str=left_join(C,lb1, by="Name")


tSplit=data.frame(matrix(NA, ncol=10, nrow=20))
colnames(tSplit)=c("seednumber","R2_test","train_RMSE","iteration",
                   "eta","maxdepth","subsample","colbytree","minchidw")
ss=sample(1:30000,60)
setwd("D:/BIOTEC NSTDA/ML_ethanol/ML_ethanol_allnew/Table")
m40c=read.csv("Mutiple_train-test_split_40C_08062023.csv")
ss=ss[!ss%in%m40c$seednumber][1:30]

for (c in 1:30) {
  print(paste("************",c,"*************"))
  
    
    ###2. prepare training and test sets ###
    set.seed(ss[c])
    A=createDataPartition(eth_str$`EtOH(g/L)`, p=0.7, list = F) #caret
    train_df=data.matrix(eth_str[A,])
    test_df=data.matrix(eth_str[-A,])
    dim(train_df);dim(test_df)
    head(test_df)
    
    ### 3. Training the model by XGBoost ###
    library(xgboost)
    
    #details of data 
    #Change number in label = n
    #set column #n 3 is ethanol,
    colnames(test_df)
    n=3
    print(colnames(test_df)[n])
    xgtrain=xgb.DMatrix(data = train_df[,c(9:11)],label=train_df[,n])
    xgtest=xgb.DMatrix(data = test_df[,c(9:11)],label=test_df[,n])
    
    
    ## Tuning ##
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
    
    # Use randomly created parameters to create 100 XGBoost-models
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
    bestrd2.1[1:5,]
    bestrd3=bestrd2.1[1,]
    bestrd3
    
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
    
    
    
    tSplit[c,1]=ss[c]
    #tSplit[c,2]=bestrd3$Rs
    tSplit[c,2]=xgb$evaluation_log[xgb$niter][,2]
    tSplit[c,3]=corr
    tSplit[c,4]=bestrd3$`mdcv$best_iteration`
    tSplit[c,5]=bestrd3$eta
    tSplit[c,6]=bestrd3$max_depth
    tSplit[c,7]=bestrd3$subsample
    tSplit[c,8]=bestrd3$colsample_bytree
    tSplit[c,9]=bestrd3$min_child_weight
    
    rm(bestrd3)
    rm(bestrd)
    rm(corr)
  }
  


setwd("D:/BIOTEC NSTDA/ML_ethanol/ML_ethanol_allnew/Table")
write.csv(tSplit, "Mutiple_train-test_split_40Cv2_09062023.csv")
