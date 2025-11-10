set.seed(12345)
library(phyloseq)
library(glmnet)
library(pROC)
library(caret)

args<-commandArgs(TRUE)
Var = "Disease"
norm = as.character(args[1]) # norm =  "Log" # scale Log
Nboot <- 1000

dir.create(norm)

#######################################################################################################################################
##########################################               Load the functions                       #####################################
#######################################################################################################################################

working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code"
 
source(paste0(working_dir,"/functions/Machine_learning_functions.R"))
source(paste0(working_dir,"/functions/supplementary_ML_functions.R"))

#######################################################################################################################################
##########################################                  Read the data                          #####################################
#######################################################################################################################################
if(norm == "Log"){
	in_phylo <- readRDS("./infiles/phylo_biomarkers_log.rds")
}
if(norm == "scale"){
	in_phylo <- readRDS("./infiles/phylo_biomarkers_scaled.rds")
}

Metadata <- data.frame(sample_data(in_phylo))
#Metadata <- Metadata[,c("Disease","enterotype","Disease_activity","BMI","Sex")]
Metadata <- Metadata[,c("Disease","enterotype","BMI","Sex","Water")]
#Metadata <- Metadata[,c("Disease","enterotype","Disease_activity","BMI","Sex","Water")]
#table(Metadata$Disease, Metadata$enterotype )
#table(Metadata$Disease_activity, Metadata$enterotype )

Metadata <- Metadata[complete.cases(Metadata),]

### Otu table ####
X_all <- data.frame(t(otu_table(in_phylo)))
X_all <- X_all[match( rownames(Metadata) , rownames(X_all) ),]

X_all$Bacteroides_2 <-   c(ifelse(Metadata$enterotype == "Bacteroides_2",1,0))
X_all$Sex <-  c(ifelse(Metadata$Sex == "Male",1,0))
X_all$y_all <-  factor(ifelse(Metadata$Disease == "SpA",1,0))

if(norm == "Log"){
	X_all$BMI <-  as.numeric(log(Metadata$BMI + 1 ))
	X_all$Water <-  as.numeric(log(Metadata$Water +1 ))
}
if(norm == "scale"){
	X_all$BMI <-  as.numeric(scale(Metadata$BMI))
	X_all$Tryptamine <-  as.numeric(scale(X_all$Tryptamine))
	X_all$Water <-  as.numeric(scale(Metadata$Water))
}
#######################################################################################################################################
##########################################             Create the outer loops                     #####################################
#######################################################################################################################################

list_cv <- list()
for(i in 1:Nboot){
	train_index <- createDataPartition(X_all$y_all, p = 0.8, list = FALSE)
	list_cv[[i]] <- train_index
	rm(train_index)
}


#######################################################################################################################################
##########################################              Run the inner loop                       ######################################
#######################################################################################################################################
importance_df <- data.frame()
stats_df <- data.frame()
list_roc <- list()
for(i in 1:length(list_cv)){
	print(i)
	train_index <- c(list_cv[[i]])
	ldf <- train_test_set(Index=train_index,X_in = X_all)

	### Model ###
	n <- nrow(ldf$X_tr)
	cv_fit <- cv.glmnet(ldf$X_tr,  ldf$Y_tr , alpha = 1, family = "binomial", nfolds = n - 1)  # o use cv.glmnet() directamente

	### Predict ####
	y_pred_prob <- predict(cv_fit, newx = as.matrix(ldf$X_val), s = "lambda.min", type = "response")
	y_pred_class <- as.numeric(ifelse(y_pred_prob > 0.5, 1, 0))

	### AUC ###
	roc_obj <- roc(ldf$Y_val, as.numeric(y_pred_prob))
	list_roc[[i]] <- roc_obj
	auc_value <- auc(roc_obj)

	conf_matrix <- confusionMatrix(factor(y_pred_class),  ldf$Y_val )
	temp_stats <- statistics_classification(confMatrix = conf_matrix)
	temp_stats$AUC <- as.numeric(auc_value)
	stats_df <- rbind(stats_df,temp_stats)
	rm(temp_stats)

	s.pecies.coef<-coef(cv_fit, s = "lambda.min")
	Coef <- as.data.frame(as.matrix(s.pecies.coef))
	ret_df <- data.frame( Varaible = Var , Predictior = rownames(Coef), Coef = Coef$lambda.min , Boot = i )
	importance_df <- rbind(  importance_df , ret_df )
	rm(ret_df,cv_fit,ldf,roc_obj)
}



#######################################
######## Plot the importance  #########
Mean = c(by(importance_df$Coef, importance_df$Predictior , mean))
SD = c(by(importance_df$Coef, importance_df$Predictior , sd))
SE = c(by(importance_df$Coef, importance_df$Predictior , standard_error))

importance_df_sig <- importance_df[importance_df$Coef != 0,]
PreFreq <- table(importance_df_sig$Predictior)


Feature = names(Mean)
MeanImp <- data.frame( Feature, Mean, SD,SE)
MeanImp$Increase <- ifelse(MeanImp$Mean>0,"Increase SpA","Decrease SpA")
MeanImp <- MeanImp[MeanImp$Feature != "(Intercept)",]
MeanImp <- MeanImp[order(MeanImp$Mean),]

MeanImp$Freq <- PreFreq[match(  as.character(MeanImp$Feature) , names(PreFreq))]
if(any(is.na(MeanImp$Freq))){
	MeanImp[is.na(MeanImp$Freq),]$Freq <- 0
}


MeanImp$Feature_Freq <- paste0( MeanImp$Feature," [",MeanImp$Freq,"/",Nboot,"]")
MeanImp$FreqCat <- ifelse(MeanImp$Freq < round(0.7 * Nboot), "Weak support", "Strong support")
MeanImp$ColorGroup <- paste(MeanImp$Increase, MeanImp$FreqCat)

MeanImp$Feature <- factor(as.character( MeanImp$Feature), as.character( MeanImp$Feature) )
MeanImp$Feature_Freq <- factor(as.character( MeanImp$Feature_Freq), as.character( MeanImp$Feature_Freq) )


my_colors <- c(  "Increase SpA Strong support" = "#B22222", "Increase SpA Weak support"  = "#F08080", 
              "Decrease SpA Strong support" = "#1E90FF", "Decrease SpA Weak support"  = "#ADD8E6" )


pImportance <- ggplot(MeanImp, aes(x = Feature_Freq, y = Mean, fill = ColorGroup)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  coord_flip() + 
  theme_bw() + 
  ylab("mean GLM Coeffcient") +
  ggtitle(paste(Var, "Feature Importance")) +
  scale_fill_manual(values = my_colors)
pImportance
ggsave(file=paste0("./",norm,"/pImportance.pdf"),pImportance, width = 10, height = 7)


mean_ROC_list <- mean_ROC_function(list_ggroc_curves = list_roc, Variable = Var)
mean_ROC_df <-  mean_ROC_list$mean_ROC 

auc_text<-paste0("AUC","\n", Var,"=",mean_ROC_list$AUC, "/", mean_ROC_list$AUCsd)

AUCMeanAUCsd <- plot_meanROC(mean_ROC_plot=mean_ROC_df, titleplot=paste("AUC",Var), auc_text_vec = auc_text, Deviation = "SD")
AUCMeanAUCsd <- AUCMeanAUCsd + labs(title = paste("AUC",Var), 
		subtitle = paste("Mean =", round(mean(stats_df$AUC,na.rm=T),digits=2) , "SD =", round(sd(stats_df$AUC,na.rm=T),digits=2) ) 
	  )	  
ggsave(file=paste0("./",norm,"/AUCMeanAUCsd.pdf"),AUCMeanAUCsd, width = 10, height = 7)


##############################
###### Save the outputs ######
save(file=paste0("./",norm,"/AUCMeanAUCsd.RData"),AUCMeanAUCsd)
save(file=paste0("./",norm,"/mean_ROC_df.RData"),mean_ROC_df)
save(file=paste0("./",norm,"/MeanImp.RData"),MeanImp)
save(file=paste0("./",norm,"/pImportance.RData"),pImportance)
save(file=paste0("./",norm,"/importance_df.RData"),importance_df)
save(file=paste0("./",norm,"/stats_df.RData"),stats_df)



