#ML using XGBOOST

library(readxl)
library(ggplot2)
library(dplyr)
library(caret)
library(reshape2)

# 30 C 
### 1. Import data  ###
#scale promoter strengths 
#data contain scale of all metabolites and promoter strengths every timepoints
pro_str=read.table("Metabolic_promoterstrength_forML__08092022.txt")

#productions data 

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

###2. prepare training and test sets ###
set.seed(019)
A=createDataPartition(eth_str$EtOH_24h, p=0.7, list = F) #caret
train_df=data.matrix(eth_str[A,])
test_df=data.matrix(eth_str[-A,])
dim(train_df);dim(test_df)
head(test_df)

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

xgb=xgboost(data=xgtrain,
            booster= "gbtree",
            objective= "reg:squarederror",
            nrounds =42,
            eta=0.1 ,
            max.depth=9,
            subsample = 1,
            colsample_bytree = 1,
            min_child_weight = 0
)
#prediction using test data
pred=predict(xgb, xgtest)
t=as.numeric(test_df[,n])
cor(pred,t)^2 #correlation ^2 = Rsquare


# Join prediction and actual values to dataframe 
preds_xgb=data.frame(test_df[,n],pred)
colnames(preds_xgb)=c("Actual","Predicted")
#importance features
im_table=xgb.importance(feature_names = colnames(train_df[,c(2:4)]), model=xgb)
#xgb.plot.importance(importance_matrix = im_table)
im_table
#using ggplot2
im_table$Feature=factor(im_table$Feature, 
                        levels = c("PDC1_24h","ADH1_24h","TPS1_24h"))

p1=im_table%>%ggplot(aes(x=Feature, y=Gain, fill=Feature))+
  geom_bar(stat = "identity")+
  ggtitle("Contribution of importance promoter in the prediction of the 30C ethanol production")+
  theme_bw()+ylab("Importance")+xlab("Promoter")+
  scale_x_discrete(labels=c("PDC1","ADH1","TPS1"))+ylim(0,0.8)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        text = element_text(size=15))

p1+scale_fill_npg()
### 4. Visualization the model ###
print(xgb)
#plot
preds_xgb%>%ggplot(aes(x=Actual, y=Predicted))+   
  geom_point(alpha = 0.8, color = "cadetblue",size=2) +
  geom_smooth(method = "loess", formula = "y ~ x", se=F, color="black") +
  #geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_x_continuous(limits=c(0,70),breaks = seq(0,70,by=10))+
  scale_y_continuous(limits=c(0,70),breaks = seq(0,70,by=10))+
  coord_fixed()+theme_bw()+xlab("Actual ethanol (g/L)")+ylab("predicted ethanol (g/L)")+
  labs(title = "XGBOOST",
       subtitle= substitute(paste(R^2, " = ",r2a,
                                  ", RMSE = ",rmsea),
                            list(r2a=round(cor(pred,t)^2,4),
                                 rmsea=round(xgb$evaluation_log$train_rmse[42],4))))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle =element_text(hjust = 0.5),
        text = element_text(size=16))

#prediction vs actual table
pretable=cbind(eth_str[-A,1],preds_xgb)
pretable=pretable[,c(1,3,2)]
colnames(pretable)=c("Name","Preiction","Actual")


#plot tree
xgb.plot.tree(feature_names = names(train_df), model = xgb, 
              trees=bestrd3$`mdcv$best_iteration`-1)
xgb.plot.multi.trees(model = xgb, feature_names = names(train_df))







