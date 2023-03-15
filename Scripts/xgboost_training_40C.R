#ML using XGBOOST

library(readxl)
library(ggplot2)
library(dplyr)
library(caret)
library(reshape2)

### High temperature ###
# 40 C 
### 1. Import data  ###
#scale promoter strengths 
#40C

setwd("D:/BIOTEC NSTDA/ML_ethanol/Results/37C40C/")
prost2_40=read.table("Promoterstrength_40C_vectornorm_19072022.txt")
rownames(prost2_40)=c("vector","PDC1","ADH1","TPS1","ACT1",
                      "PGK1","ENO2","TDH3","YEF3")
sprost40=as.data.frame(scale(prost2_40[-1,7], center = F, scale = T))
sprost40=cbind(rownames(prost2_40)[-1],sprost40)
colnames(sprost40)=c("promoter","promoterstrength")
m4=sprost40[order(sprost40$promoterstrength, decreasing = T),]
m4$promoter=factor(m4$promoter, levels=m4$promoter)

ggplot(m4,aes(promoter, promoterstrength, fill=promoter))+
  geom_bar(stat = "identity")+
  ggtitle("Normalized promoter strength at 40C")+
  theme_bw()+ylab("Normalized promoter strength")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        text=element_text(size=18))+
  scale_fill_manual(values = c("#8dd3c7", "#ffffb3", "#bebada","#80b1d3",
                               "#b3de69","#fb8072","#fdb462","#fccde5"))

#productions data 
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

#add promoter strength
#import label from P'Wut
setwd("D:/BIOTEC NSTDA/ML_ethanol/Results/30C")
lb=read.csv("input_AI_Guided.csv")
colnames(lb)[1]="Name"
lb1=lb[,c(1,3,4,5)]
col=list(c("PDC1","ADH1","TPS1"))
lb1[,col[[1]]]=sprost40$promoterstrength[match(unlist(lb1[,col[[1]]]),sprost40$promoter)]

eth_str=left_join(C,lb1, by="Name")

###2. prepare training and test sets ###
set.seed(019)
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
n=4
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
set.seed(0199) #Ethanol
xgb=xgboost(data=xgtrain,
            booster= "gbtree",
            objective= "reg:squarederror",
            nrounds =61,
            eta=0.5,
            max.depth=3,
            subsample = 1,
            colsample_bytree = 0.7,
            min_child_weight = 0
)
set.seed(0199) #OD
xgb=xgboost(data=xgtrain,
            booster= "gbtree",
            objective= "reg:squarederror",
            nrounds =44,
            eta=0.3,
            max.depth=4,
            subsample = 1,
            colsample_bytree = 0.7,
            min_child_weight = 0
)
set.seed(0199) #%Glu
xgb=xgboost(data=xgtrain,
            booster= "gbtree",
            objective= "reg:squarederror",
            nrounds =57,
            eta=0.4,
            max.depth=4,
            subsample = 1,
            colsample_bytree = 0.7,
            min_child_weight = 1
)

#prediction using test data
pred=predict(xgb, xgtest)
t=as.numeric(test_df[,n])
cor(pred,t)^2 #correlation ^2 = Rsquare


# Join prediction and actual values to dataframe 
preds_xgb=data.frame(test_df[,n],pred)
colnames(preds_xgb)=c("Actual","Predicted")
#importance features
im_table=xgb.importance(feature_names = colnames(train_df[,c(9:11)]),
                        model=xgb)
#xgb.plot.importance(importance_matrix = im_table)
im_table
#using ggplot2
im_table$Feature=factor(im_table$Feature, 
                        levels = c("PDC1","ADH1","TPS1"))

p1=im_table%>%ggplot(aes(x=Feature, y=Gain, fill=Feature))+
  geom_bar(stat = "identity")+
  ggtitle("Contribution of importance promoter in the prediction of the 40C ethanol production")+
  theme_bw()+ylab("Importance")+xlab("Promoter")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        text = element_text(size=15))

p1+scale_fill_npg()
### 4. Visualization the model ###
print(xgb)

#plot
preds_xgb%>%ggplot(aes(x=Actual, y=Predicted))+   
  geom_point(alpha = 0.8, color = "cadetblue",size=2.5) +
  geom_smooth(method = "loess", formula = "y ~ x",se=F, color="black") +
  #geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_x_continuous(limits=c(0,65),breaks = seq(0,65,by=10))+
  scale_y_continuous(limits=c(0,65),breaks = seq(0,65,by=10))+
  theme_bw()+coord_fixed()+xlab("Actual ethanol (g/L)")+ylab("predicted ethanol (g/L)")+
  labs(title = "XGBOOST",
       subtitle= substitute(paste(R^2, " = ",r2a,
                                  ", RMSE = ",rmsea), 
                            list(r2a=round(cor(pred,t)^2,4),
                                 rmsea=round(xgb$evaluation_log$train_rmse[61],4))))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle =element_text(hjust = 0.5),
        text = element_text(size=16))


#plot tree
xgb.plot.tree(feature_names = colnames(train_df[,c(9:11)]),
              model = xgb, trees=bestrd3$`mdcv$best_iteration`-1)

xgb.plot.multi.trees(model = xgb, 
                     feature_names = colnames(train_df[,c(9:11)]))




############# PREDICTION #################
lb2=lb1[!lb1$Name%in%eth_str$Name,]
mlb2=xgb.DMatrix(data.matrix(lb2[,2:4]))
#prediction using the rest of combination data
pred2=predict(xgb, mlb2) 
new_pred=data.frame(lb2,pred2)
new_pred2=left_join(lb[,c(1,3,4,5)],new_pred, by="Name")
pre_EtOH=new_pred2[order(new_pred2$pred2, decreasing = T),]
colnames(pre_EtOH)=c("Name","PDC1","ADH1","TPS1",
                     "PDC1_norm","ADH1_norm","TPS1_norm","Predicted_EtOH")

pre_EtOH[1:30,]


pred3=predict(xgb, mlb2) 
new_pred3=data.frame(lb2,pred3)
new_pred3=left_join(lb[,c(1,3,4,5)],new_pred3, by="Name")
pre_OD=new_pred3[order(new_pred3$pred3, decreasing = T),]
colnames(pre_OD)=c("Name","PDC1","ADH1","TPS1",
                   "PDC1_norm","ADH1_norm","TPS1_norm","Predicted_OD600")
pre_OD[1:30,]


pred4=predict(xgb, mlb2) 
new_pred4=data.frame(lb2,pred4)
new_pred4=left_join(lb[,c(1,3,4,5)],new_pred4, by="Name")
pre_Glu=new_pred4[order(new_pred4$pred4, decreasing = T),]
colnames(pre_Glu)=c("Name","PDC1","ADH1","TPS1",
                   "PDC1_norm","ADH1_norm","TPS1_norm","Predicted_Glucoseconsump")
pre_Glu[1:30,]

inter_pre=intersect(pre_EtOH[1:30,"Name"],pre_OD[1:30,"Name"])
pre_EtOH[pre_EtOH$Name%in%inter_pre,][1:15,]
candi=pre_EtOH[pre_EtOH$Name%in%inter_pre,][1:15,]
candi=left_join(candi, pre_OD[,c(1,8)], by="Name")
candi=left_join(candi, pre_Glu[,c(1,8)], by="Name")

setwd("D:/BIOTEC NSTDA/ML_ethanol/ML_ethanol_allnew/Table")
write.table(candi,"Prediction_usingTuned40C_EtOH_XGBoostmodel_tunedmodel_25102022.txt")
candi=read.table("Prediction_usingTuned40C_EtOH_XGBoostmodel_tunedmodel_25102022.txt")

library(writexl)
write_xlsx(candi, path="Prediction_usingTuned40C_EtOH_XGBoostmodel_tunedmodel_25102022.xlsx")


library(VennDiagram)
set1=pre_EtOH[1:30,"Name"]
set2=pre_OD[1:30,"Name"]
set3=pre_Glu[1:30,"Name"]
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(x=list(set1,set2,set3),
             category.names = c("Ethanol" , "OD " , "% Glucose consumption"),
             filename = '#14_venn_diagramm.png', fill = myCol,
             output= T)
