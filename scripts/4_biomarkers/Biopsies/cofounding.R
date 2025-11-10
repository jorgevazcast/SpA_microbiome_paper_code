set.seed(12345)
library(glmmTMB)    
library(glmnet)
library(lmerTest)
library(microbiome)

argscomd = commandArgs(trailingOnly=TRUE)
print(argscomd[1])

GMPR.file <- as.character(argscomd[1])
rar.file <- as.character(argscomd[2])
above10000.file <- as.character(argscomd[3])
clr.file <- as.character(argscomd[4])
Variable = as.character(argscomd[5]) ## Variable = "Disease_activity" # Histological_Inflammation Histology Disease
Tax = as.character(argscomd[6]) ## "GreenGenes2" "GTDB_r86" "GTDB_r202" "GTDB_r207" "GTDB_r220" "rdp_16" "rdp_19" "silva_v138_2" "Consensus_sp"
sig_features =  as.character(argscomd[7]) # 

cat("\n Parameters: ", "\n", GMPR.file, "\n", rar.file, "\n", above10000.file, "\n", clr.file,"\n", Variable,"\n", Tax,"\n",sig_features,"\n" )

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

### Read the taxonomy
Sp_annotations <- readRDS(paste0(working_dir,"/1_infiles/Colon_Ileum_Biopsy/Sp_annotations.rds"))

# Taxonomy;: "GreenGenes2" "GTDB_r86" "GTDB_r202" "GTDB_r207" "GTDB_r220" "rdp_16" "rdp_19" "silva_v138_2" "Consensus_sp"
Annotation <- ASV_naming_function( Sp_annotations=Sp_annotations, Taxonomy=Tax)

### Read the phyloseq objects
physeq.GMPR <- read_input_phyloseq(GMPR.file)
physeq.rar <- read_input_phyloseq(rar.file)
physeq.clr <- read_input_phyloseq(clr.file)
physeq.all <- read_input_phyloseq(above10000.file)

### Update the taxonomy
physeq.GMPR <- phyloseq_asv_annotation(in.phylo = physeq.GMPR , Annotation = Annotation)
physeq.rar <- phyloseq_asv_annotation(in.phylo = physeq.rar , Annotation = Annotation)
physeq.clr <- phyloseq_asv_annotation(in.phylo = physeq.clr , Annotation = Annotation)
physeq.all <- phyloseq_asv_annotation(in.phylo = physeq.all , Annotation = Annotation)

#######################################################################################################################################################
##########################################                  Read the significant features                        ######################################
#######################################################################################################################################################
q_values_table <- read.table(sig_features,sep="\t", header = T)

q_values_table <- subset(q_values_table, ACAT <= p_val_cut_off ) ### The ACAT must be significant
Sig_taxas <- apply(   q_values_table[, 2:(ncol(q_values_table)-1)], 1, function(x) { any(x <= p_val_cut_off)  } ) ### Other test must be significat
q_values_table_sig <- q_values_table[Sig_taxas,]


############################################################################################################################################
##########################################                   Cofound by the age                   ##########################################
############################################################################################################################################

Features <- unique(q_values_table_sig$Feature)

clr_cofound_test_anova <- data.frame()
clr_cofound_summary <- data.frame()
for(i in  Features){

	subphylo_df <- psmelt( prune_taxa(i,physeq.clr) )
	tempList <- glm_cofound_function(Y.var = "Abundance", NULL_model = c("Age_at_visit"), Test_model = c("Age_at_visit",Variable), in.df = subphylo_df)  

	stast_anova <- data.frame(Feature = i, tempList[["Stats"]])
	clr_cofound_test_anova <- rbind(clr_cofound_test_anova,stast_anova)
	
	temp_summary <- summary(tempList[["Model"]])
	temp_summary <- data.frame(temp_summary$coefficients)
	temp_summary <- data.frame(Feature = i, Variable = rownames(temp_summary) ,temp_summary)
	clr_cofound_summary <- rbind(clr_cofound_summary, temp_summary)

	rm(stast_anova,temp_summary,subphylo_df)
}
rownames(clr_cofound_summary) <- NULL
clr_cofound_test_anova
clr_cofound_test_anova$q.value <- p.adjust(clr_cofound_test_anova$Pr_Chi,method="BH")
clr_cofound_summary$q.value  <- p.adjust(clr_cofound_summary$Pr...t..,method="BH")

write.table(file="clr_cofound_test_anova.tsv",clr_cofound_test_anova,sep = "\t", col.names = TRUE, row.names = F)
write.table(file="clr_cofound_summary.tsv",clr_cofound_summary,sep = "\t", col.names = TRUE, row.names = F)


gmpr_cofound_test_anova <- data.frame()
gmpr_cofound_summary <- data.frame()
for(i in  Features){

	subphylo_df <- psmelt( prune_taxa(i,physeq.GMPR) )
	tempList <- nb_cofound_function(Y.var = "Abundance", NULL_model = c("Age_at_visit"), Test_model = c("Age_at_visit",Variable), in.df = subphylo_df)  

	stast_anova <- data.frame(Feature = i, tempList[["Stats"]])
	gmpr_cofound_test_anova <- rbind(gmpr_cofound_test_anova,stast_anova)
	
	temp_summary <- summary(tempList[["Model"]])
	temp_summary <- data.frame(temp_summary$coefficients$cond)
	temp_summary <- data.frame(Feature = i, Variable = rownames(temp_summary) ,temp_summary)
	gmpr_cofound_summary <- rbind(gmpr_cofound_summary, temp_summary)

	rm(stast_anova,temp_summary,subphylo_df)
}
rownames(gmpr_cofound_summary) <- NULL
gmpr_cofound_test_anova
gmpr_cofound_test_anova$q.value <- p.adjust(gmpr_cofound_test_anova$Pr_Chi,method="BH")
gmpr_cofound_summary$q.value  <- p.adjust(gmpr_cofound_summary$Pr...z..,method="BH")


write.table(file="gmpr_cofound_summary.tsv",gmpr_cofound_summary,sep = "\t", col.names = TRUE, row.names = F)
write.table(file="gmpr_cofound_test_anova.tsv",gmpr_cofound_test_anova,sep = "\t", col.names = TRUE, row.names = F)



