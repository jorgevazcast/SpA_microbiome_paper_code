set.seed(12345)
library(glmmTMB)    
library(glmnet)
library(lmerTest)
library(microbiome)

argscomd = commandArgs(trailingOnly=TRUE)

log.file <- as.character(argscomd[1])
Variable = as.character(argscomd[2]) 
sig_features =  as.character(argscomd[3])

cat("\n Parameters: ", "\n", log.file,"\n", Variable,"\n", sig_features,"\n" )

###############################################
######## Create the outfiles directory ########

dir2create <- paste0(getwd(),"/",Variable)
dir.create(file.path(dir2create ), showWarnings = FALSE)
setwd(file.path(dir2create))

#################################
######## q-value cut-off ########
p_val_cut_off <- 0.1

#######################################################################################################################################
##########################################              Load the functions                       ######################################
#######################################################################################################################################

working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

source(paste0(working_dir,"/functions/data_processing_functions.R"))
source(paste0(working_dir,"/functions/statistical_test_functions.R"))
source(paste0(working_dir,"/functions/CATDB_supplementary_functions.R"))
source(paste0(working_dir,"/functions/Statistical_model_tools.R"))

#######################################################################################################################################
##########################################                  Read the data                        ######################################
#######################################################################################################################################

### Read the phyloseq objects
physeq.log <- read_input_phyloseq(log.file)

taxa_names(physeq.log) <- gsub("[^[:alnum:]]", "_", taxa_names(physeq.log))
taxa_names(physeq.log) <- gsub(" ", "_", taxa_names(physeq.log))	

#######################################################################################################################################################
##########################################                  Read the significant features                        ######################################
#######################################################################################################################################################

q_values_table <- read.table(sig_features,sep="\t", header = T)

q_values_table <- subset(q_values_table, ACAT <= p_val_cut_off ) ### The ACAT must be significant
Sig_taxas <- apply(   q_values_table[, 2:(ncol(q_values_table)-1)], 1, function(x) { any(x <= p_val_cut_off)  } ) ### Other test must be significat
q_values_table_sig <- q_values_table[Sig_taxas,]

#######################################################################################################################################################
##########################################               Cofound by the significant features                     ######################################
#######################################################################################################################################################

Features <- unique(q_values_table_sig$Feature)
log_cofound_test_anova <- data.frame()
log_cofound_summary <- data.frame()
for(i in  Features){

	subphylo_df <- psmelt( prune_taxa(i,physeq.log) )
	tempList <- glm_cofound_function(Y.var = "Abundance", NULL_model = c("Water","Age_at_visit"), 
			Test_model = c("Water","Age_at_visit",Variable), in.df = subphylo_df)  

	stast_anova <- data.frame(Feature = i, tempList[["Stats"]])
	log_cofound_test_anova <- rbind(log_cofound_test_anova,stast_anova)
	
	temp_summary <- summary(tempList[["Model"]])
	temp_summary <- data.frame(temp_summary$coefficients)
	temp_summary <- data.frame(Feature = i, Variable = rownames(temp_summary) ,temp_summary)
	log_cofound_summary <- rbind(log_cofound_summary, temp_summary)

	rm(stast_anova,temp_summary,subphylo_df)
}
rownames(log_cofound_summary) <- NULL
#log_cofound_test_anova
log_cofound_test_anova$q.value <- p.adjust(log_cofound_test_anova$Pr_Chi,method="BH")
log_cofound_summary$q.value  <- p.adjust(log_cofound_summary$Pr...t..,method="BH")


write.table(file="log_cofound_summary.tsv",log_cofound_summary,sep = "\t", col.names = TRUE, row.names = F)
write.table(file="log_cofound_test_anova.tsv",log_cofound_test_anova,sep = "\t", col.names = TRUE, row.names = F)



