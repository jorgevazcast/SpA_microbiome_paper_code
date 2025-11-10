set.seed(12345)
library(phyloseq)
library(glmnet)
library(pROC)
library(caret)

args<-commandArgs(TRUE)
Var = "Response"
norm = as.character(args[1]) # norm =  "Log" # scale Log none
Nboot <- 1000
q_val_cut_off <- 0.1

dir.create(norm)

working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

biomarkersFile <- paste0(working_dir,"/4_biomarkers/mOTUs_gogut/DA/Response/q_values_table.tsv")

#######################################################################################################################################
##########################################               Load the functions                       #####################################
#######################################################################################################################################

source(paste0(working_dir,"/functions/Machine_learning_functions.R"))
source(paste0(path_project,"/functions/supplementary_ML_functions.R"))

#######################################################################################################################################
##########################################                  Read the data                          ####################################
#######################################################################################################################################
if(norm == "Log"){
	in_phylo <- readRDS("./infiles/phylo_biomarkers_log.rds")
}else if(norm == "scale"){
	in_phylo <- readRDS("./infiles/phylo_biomarkers_scaled.rds")
}else if(norm == "none"){
	in_phylo <- readRDS("./infiles/phylo_biomarkers.rds")
}

q_values_table <- read.table(biomarkersFile,header = T,sep="\t")
q_values_table <- subset(q_values_table, ACAT <= q_val_cut_off ) ### The ACAT must be significant
Sig_taxas <- apply(   q_values_table[, 2:(ncol(q_values_table)-1)], 1, function(x) { any(x <= q_val_cut_off)  } ) ### Other test must be significat
q_values_table_sig <- q_values_table[Sig_taxas,]
in_phylo <- prune_taxa( q_values_table_sig$Feature, in_phylo)

Metadata <- data.frame(sample_data(in_phylo))
Metadata <- Metadata[complete.cases(Metadata$Response),]

### Otu table ####
X_all <- data.frame(t(otu_table(in_phylo)))
X_all <- X_all[match( rownames(Metadata) , rownames(X_all) ),]

X_all$y_all <-  factor(ifelse(Metadata$Response == "Responder",1,0) )

#######################################################################################################################################
##########################################              Run the inner loop                       ######################################
#######################################################################################################################################
importance_df <- data.frame()
sum <- 0
while (sum <= Nboot) {
#for(i in 1:Nboot){

	
	ldf <- sample( 1:nrow(X_all) , size = nrow(X_all), replace = T )
	Xboot <- X_all[ldf,]
	Y_tr <- Xboot$y_all	
	Xboot$y_all <- NULL
	if( table(Y_tr)["0"] > 2 ){
	 
		print(table(Y_tr)["0"])
		### Model ###
		n <- nrow(Xboot)
		cv_fit <- cv.glmnet(as.matrix(Xboot),  Y_tr , alpha = 1, family = "binomial", nfolds = c(n - 1) )  # 

		### Predict ####
		s.pecies.coef<-coef(cv_fit, s = "lambda.min")
		Coef <- as.data.frame(as.matrix(s.pecies.coef))
		ret_df <- data.frame( Varaible = Var , Predictior = rownames(Coef), Coef = Coef$lambda.min , Boot = sum )
		importance_df <- rbind(  importance_df , ret_df )
		sum <- sum + 1
		print(sum)
		
	}
	rm(ret_df,cv_fit,ldf)
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
MeanImp$Increase <- ifelse(MeanImp$Mean >0 ,"Responder","Non_responder")
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


my_colors <- c(  "Non_responder Strong support" = "#B22222", "Non_responder Weak support"  = "#F08080", 
              "Responder Strong support" = "#1E90FF", "Responder Weak support"  = "#ADD8E6" )

pImportance <- ggplot(MeanImp, aes(x = Feature_Freq, y = Mean, fill = ColorGroup)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  coord_flip() + 
  theme_bw() + 
  ylab("mean GLM Coeffcient") +
  ggtitle(paste(Var, "Feature Importance")) +
  scale_fill_manual(values = my_colors)
pImportance
ggsave(file=paste0("./",norm,"/pImportance.pdf"),pImportance, width = 12, height = 6)

##############################
###### Save the outputs ######
save(file=paste0("./",norm,"/MeanImp.RData"),MeanImp)
save(file=paste0("./",norm,"/pImportance.RData"),pImportance)
save(file=paste0("./",norm,"/importance_df.RData"),importance_df)

