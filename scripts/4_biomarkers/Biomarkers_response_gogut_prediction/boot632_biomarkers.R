set.seed(12345)
library(phyloseq)
library(glmnet)
library(pROC)
library(caret)
library(DMwR)
library(fastshap)
library(shapviz)

args<-commandArgs(TRUE) # 
Var = as.character(args[1]) #   Var = "Response"   
norm = as.character(args[2]) # norm = "Log" # norm = "scale" Log  norm = "none"
method = as.character(args[3]) # method = "glmnet" method = "rf"
include_biomarkers = as.logical(as.character(args[4])) #  include_biomarkers = T   include_biomarkers =  F
Nboot <- 1000
q_val_cut_off <- 0.1

working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

biomarkersFile <- paste0(working_dir,"/4_biomarkers/mOTUs_gogut/DA/Response/q_values_table.tsv")

outdir <- paste0("boot632_",method,"_",norm,"_","only_FS_",include_biomarkers,"_",Var)
dir.create(outdir)

#######################################################################################################################################
##########################################               Load the functions                       #####################################
#######################################################################################################################################

source(paste0(working_dir,"/Functions/Machine_learning_functions.R"))
source(paste0(working_dir,"/functions/supplementary_ML_functions.R"))

#######################################################################################################################################
##########################################                  Read the data                          #####################################
#######################################################################################################################################
if(norm == "Log"){
	in_phylo <- readRDS("./infiles/phylo_biomarkers_log.rds")
} else if(norm == "scale"){
	in_phylo <- readRDS("./infiles/phylo_biomarkers_scaled.rds")
} else if(norm == "none"){
	in_phylo <- readRDS("./infiles/phylo_biomarkers.rds")
}



if(include_biomarkers == T){
	q_values_table <- read.table(biomarkersFile,header = T,sep="\t")
	q_values_table <- subset(q_values_table, ACAT <= q_val_cut_off ) ### The ACAT must be significant
	Sig_taxas <- apply(   q_values_table[, 2:(ncol(q_values_table)-1)], 1, function(x) { any(x <= q_val_cut_off)  } ) ### Other test must be significat
	q_values_table_sig <- q_values_table[Sig_taxas,]
	in_phylo <- prune_taxa( q_values_table_sig$Feature, in_phylo)
}else{
	PrevFilt <- prevalence(in_phylo)
	PrevFilt <- names(PrevFilt[PrevFilt >= 0.2])
	in_phylo <- prune_taxa( PrevFilt, in_phylo)	
}

Metadata <- data.frame(sample_data(in_phylo))
Metadata$ASDAS_W0 <- as.numeric(as.character(Metadata$ASDAS_W0))
Metadata$Fecal_calpro_values <- as.numeric(as.character(Metadata$Fecal_calpro_values))

Metadata$Var <- factor(as.character(Metadata$Response) ,c("Responder","Non_responder") )
Metadata <- Metadata[,c("Var","Bacteroides_2","Fecal_calpro_values","ASDAS_W0","Bacteroides_1","Ruminococcus","Prevotella" )]
Metadata <- Metadata[complete.cases(Metadata),]

### Otu table ####
X_all <- data.frame(t(otu_table(in_phylo)))
X_all <- X_all[match( rownames(Metadata) , rownames(X_all) ),]

X_all$Bacteroides_1 <-   Metadata$Bacteroides_1
X_all$Bacteroides_2 <-   Metadata$Bacteroides_2
X_all$Prevotella <-   Metadata$Prevotella
X_all$Ruminococcus <-   Metadata$Ruminococcus

X_all$y_all <-  Metadata$Var

if(norm == "Log"){
	X_all$Fecal_calpro_values <-  as.numeric(log( Metadata$Fecal_calpro_values +1 ))
	X_all$ASDAS_W0 <-  as.numeric(log( Metadata$ASDAS_W0 +1 ))
		
}else if(norm == "scale"){
	X_all$Fecal_calpro_values <-  as.numeric(scale(Metadata$Fecal_calpro_values))
	X_all$ASDAS_W0 <-  as.numeric(scale(Metadata$ASDAS_W0))

}else if(norm == "none"){
	X_all$Fecal_calpro_values <- as.numeric( Metadata$Fecal_calpro_values )
	X_all$ASDAS_W0 <- as.numeric( Metadata$ASDAS_W0 )

}
#######################################################################################################################################
##########################################             Parameters for the ml run                     ###################################
#######################################################################################################################################

# Crear trainControl
ctrl <- trainControl(
  method = "boot632",     # Bootstrap .632
  number = 1000,           #  N bootstraps
#  sampling = "smote",     # SMOTE per iteration
  classProbs = TRUE,      # for the AUC estimation
  summaryFunction = twoClassSummary,  # For the estimation of the AUC, Sens, Spec, etc.
  savePredictions = "final",  # safe the predictons
  verboseIter = F
)

#######################################################################################################################################
##########################################              Run the outer loop                       ######################################
#######################################################################################################################################

X <- X_all[, !(colnames(X_all) %in% "y_all")]
y <- X_all$y_all

#  glmnet (Lasso/Elastic Net) or RF
fit <- train(
	x = X,
	y = y,
	method = method,
	metric = "ROC",
	trControl = ctrl,
	tuneLength = 10  
)

# Resultados
print(fit)
plot(fit)

saveRDS(file=paste0("./",outdir,"/fit.rds"),fit)

#######################################################################################################################################
##########################################              Plot the results                        ######################################
#######################################################################################################################################

# 1. Extract OOB predictions from all replicates
oob_predictions <- fit$pred  # Contains columns: 'obs', 'case', 'control', 'rowIndex', etc.

best_hyper <- fit$bestTune

if(method == "glmnet"){
# 2. Filter predictions for the final model (best hyperparameters)

	final_predictions <- subset(oob_predictions,  alpha == best_hyper$alpha &  lambda == best_hyper$lambda)
} else if(method == "rf"){
# For random forest models you would use mtry instead:
	final_predictions <- subset(oob_predictions, mtry == best_hyper$mtry)
}

# 3. Calculate and plot ROC
roc_curve <- roc(response = final_predictions$obs,
                 predictor = final_predictions$Responder,  # Use probability for "Responder" class
                 levels = c("Responder", "Non_responder"))  # Note level ordering!

pdf(paste0("./",outdir,"/ROC_curve.pdf"),width = 7, height = 6)
plot(roc_curve, main = "ROC Curve (Bootstrap .632)", col = "blue", print.auc = TRUE,  print.auc.x = 0.2,  print.auc.y = 0.2)
dev.off()
auc(roc_curve)  # Print AUC value

# Calculate variable importance
importance <- varImp(fit, scale = FALSE)

# Convert to dataframe and clean names
importance_df <- as.data.frame(importance$importance)
importance_df$Variable <- rownames(importance_df)
rownames(importance_df) <- NULL
colnames(importance_df)[1] <- "Importance"

# Order variables by importance
importance_df <- importance_df[order(-importance_df$Importance), ]
importance_df<- importance_df[importance_df$Importance != 0,]

# Create the plot
plotImp <- ggplot(importance_df, aes(x = reorder(Variable, Importance), 
		y = Importance, 
		fill = Importance)) +
		geom_col(width = 0.7) +
		scale_fill_gradient(low = "steelblue1", high = "steelblue4") +
		coord_flip() +
		labs(title = "Unscaled Variable Importance",
			subtitle = paste("Model:", fit$method),
			x = "Predictors",
			y = "Importance (Absolute Coefficient)") +
		theme_classic() +
		theme(legend.position = "none",
			plot.title = element_text(face = "bold", size = 14),
			axis.text.y = element_text(size = 10),
		panel.grid.major.y = element_blank()) +
		geom_text(aes(label = round(Importance, 2)), 
		hjust = -0.1, size = 3.5)


if(include_biomarkers == T){
	ggsave(file = paste0("./",outdir,"/plotImp.pdf"), plotImp , width = 12, height = 6 )
}else if(include_biomarkers == F){
	ggsave(file = paste0("./",outdir,"/plotImp.pdf"), plotImp , width = 12, height = 20 )

}

#	ggsave(file = paste0("./",outdir,"/plotImp.pdf"), plotImp , width = 12, height = 5 )

###############################################
####### Obtain the model coeffcients ##########

final_glmnet <- fit$finalModel  

# Optimal lambda
best_lambda <- fit$bestTune$lambda
best_lambda

# Extraer coeficientes para ese lambda
coef_matrix <- coef(final_glmnet, s = best_lambda)

# Pasar a un data.frame con las features no nulas
selected_features <- data.frame(
  Feature = rownames(coef_matrix),
  Coefficient = as.numeric(coef_matrix)
)

# Filtrar solo las especies seleccionadas (coef ≠ 0)
selected_features <- selected_features[selected_features$Coefficient != 0, ]
selected_features <- selected_features[order(selected_features$Coefficient),]
selected_features$Variable <- factor(selected_features$Feature , as.character(selected_features$Feature) )

selected_features$Response <- ifelse( selected_features$Coefficient > 0 , "Non_responder","Responder" )

plot_glmnet <- ggplot(selected_features, aes(x = Variable, y = Coefficient, fill = Response)) +
	geom_col(width = 0.7) +
	coord_flip() +
	labs(title = "glmnet coefficients",
		subtitle = paste("Model:", fit$method),
		x = "Predictors",
		y = "Coefficient") +
	theme_classic() + 
	theme(plot.title = element_text(face = "bold", size = 14),axis.text.y = element_text(size = 10),panel.grid.major.y = element_blank()) +
	geom_text(aes(label = round(Coefficient, 2)), hjust = -0.1, size = 3.5) + 
	scale_fill_manual(values = c("red", "blue" ))
plot_glmnet

if( dim(selected_features)[1] < 20){
	ggsave(file = paste0("./",outdir,"/glmnet_coef.pdf"), plot_glmnet, width = 14, height = 7)
} else if( dim(selected_features)[1] > 20){
	ggsave(file = paste0("./",outdir,"/glmnet_coef.pdf"), plot_glmnet, width = 14, height = 10)
} else if( dim(selected_features)[1] > 100){
	ggsave(file = paste0("./",outdir,"/glmnet_coef.pdf"), plot_glmnet, width = 14, height = 17)
} else if( dim(selected_features)[1] > 160){
	ggsave(file = paste0("./",outdir,"/glmnet_coef.pdf"), plot_glmnet, width = 14, height = 20)
}

#######################################################################################################################################
##########################################                 Run fastshap                          ######################################
#######################################################################################################################################

###############################
####### Run fastshap ##########

# Get the best lambda chosen by caret
best_lambda <- fit$bestTune$lambda  

# Prediction wrapper for glmnet
pred_fun <- function(object, newdata, ...) {
	# Predict probabilities (type = "response") using the best lambda
	as.numeric(
		predict(
			object,
			newx = as.matrix(newdata),
			s = best_lambda,
			type = "response"
		)
	)
}

# Compute SHAP values
shap_values <- fastshap::explain(
	object = fit$finalModel,	# glmnet model
	X = X,						# feature matrix
	pred_wrapper = pred_fun,	# custom prediction function
	nsim = 1000					# number of Monte Carlo simulations
)

# Preview SHAP values
head(shap_values)

# 4. Create shapviz object
# X = human-readable feature dataframe
sv <- shapviz(shap_values, X = as.data.frame(X))

# Global importance (mean |SHAP|)
Barplot <- sv_importance(sv, kind = "bar")

# Beeswarm plot (shows direction and distribution)
sv_importance(sv, kind = "beeswarm")

# Dependence plot for one feature
sv_dependence(sv, v = names(X)[1])

saveRDS(file=paste0("./",outdir,"/SHAP_object.rds"),sv)

###############################
####### Plot fastshap ##########

Join_data <- Barplot$data
Join_data$Response <- selected_features[match( rownames(Join_data) , selected_features$Feature),]$Response
plot_SHAP_glmnet <- ggplot(Join_data, aes(x = feature, y = value, fill = Response)) +
	geom_col(width = 0.7) +
	coord_flip() +
	labs(title = "SHAP value",
		subtitle = paste("Model:", fit$method),
		x = "Predictors",
		y = "mean(|SHAP value|)") +
	theme_classic() + 
	theme(plot.title = element_text(face = "bold", size = 14),axis.text.y = element_text(size = 10),panel.grid.major.y = element_blank()) +
	geom_text(aes(label = round(value, 2)), hjust = -0.1, size = 3.5) + 
	scale_fill_manual(values = c("red", "blue" ))
plot_glmnet

if( dim(Barplot$data)[1] < 20){
	ggsave(file = paste0("./",outdir,"/SHAP_Imp.pdf"), Barplot, width = 10, height = 7)
	ggsave(file = paste0("./",outdir,"/SHAP_glmnet.pdf"), plot_SHAP_glmnet, width = 12, height = 7)
} else if( dim(Barplot$data)[1] > 20){
	ggsave(file = paste0("./",outdir,"/SHAP_Imp.pdf"), Barplot, width = 10, height = 10)
	ggsave(file = paste0("./",outdir,"/SHAP_glmnet.pdf"), plot_SHAP_glmnet, width = 12, height = 10)	
} else if( dim(Barplot$data)[1] > 100){
	ggsave(file = paste0("./",outdir,"/SHAP_Imp.pdf"), Barplot, width = 10, height = 17)
	ggsave(file = paste0("./",outdir,"/SHAP_glmnet.pdf"), plot_SHAP_glmnet, width = 12, height = 17)	
} else if( dim(Barplot$data)[1] > 160){
	ggsave(file = paste0("./",outdir,"/SHAP_Imp.pdf"), Barplot, width = 10, height = 20)
	ggsave(file = paste0("./",outdir,"/SHAP_glmnet.pdf"), plot_SHAP_glmnet, width = 12, height = 20)	
}

saveRDS(file=paste0("./",outdir,"/SHAP_glmnet_data_frame.rds"),Join_data)


