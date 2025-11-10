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


#######################################################################################################################################
##########################################                  Read the data                        ######################################
#######################################################################################################################################

### Read the taxonomy
Sp_annotations <- readRDS(paste0(working_dir,"/1_infiles/Colon_Ileum_Biopsy/Sp_annotations.rds"))

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


#######################################################################################################################################################
##########################################                   Cofound by the age                         ######################################
#######################################################################################################################################################

Variables <- c("Disease", "Diagnosis", "Disease_activity", "NSAID_use", "Age_at_visit", "Sex", "Colon_inflammation", "Ileum_inflammation", "gut_inflammation", "artritis", "uveitis", "dactylitis", "inflammatory_back_pain", "enthesitis_general", "psoriasis_skin")

Features <- unique(q_values_table_sig$Feature)

###################################################
#####################   CLR   #####################

#### Test the individual features using the model and usign the Anova function to check if its significant
clr_test_anova <- data.frame()
for(i in  Features){

	subphylo_df <- psmelt( prune_taxa(i,physeq.clr) )	
	temp_res <- data.frame()
	for(j in Variables){
		tempList <- glm_function(Y.var = "Abundance", X.var = j, in.df=subphylo_df )
		tempAnova <- data.frame(Feature = i,tempList$anova)
		temp_res <- rbind(temp_res,tempAnova)
		rm(tempList,tempAnova)
	}
	temp_res$q.value <- p.adjust(temp_res$Pr..Chisq., method = "BH")
	clr_test_anova <- rbind(clr_test_anova,temp_res)
	rm(temp_res,subphylo_df)
}

clr_test_anova <- clr_test_anova[order(clr_test_anova$Pr..Chisq.),]

write.table(file="clr_test_anova.tsv",clr_test_anova,sep = "\t", col.names = TRUE, row.names = F)
sig_clr_test_anova <- clr_test_anova[clr_test_anova$q.value <= 0.1,]

#### Cofound the analysis ####
clr_cofound_test_anova <- data.frame()
clr_cofound_summary <- data.frame()
for(i in unique(sig_clr_test_anova$Feature)){

	tempsig <- sig_clr_test_anova[sig_clr_test_anova$Feature == i,]
	subphylo_df <- psmelt( prune_taxa(i,physeq.clr) )	
	
	for(j in unique(tempsig$var)){
		tempList <- glm_cofound_function(Y.var = "Abundance", NULL_model = c("Age_at_visit"), 
				Test_model = c("Age_at_visit",j), in.df = subphylo_df)  
		stast_anova <- data.frame(Feature = i,Variable2Test = j, tempList[["Stats"]])
		clr_cofound_test_anova <- rbind(clr_cofound_test_anova,stast_anova)
	
		temp_summary <- summary(tempList[["Model"]])
		temp_summary <- data.frame(temp_summary$coefficients)
		temp_summary <- data.frame(Feature = i, Variable2Test = j, Varaible = rownames(temp_summary) ,temp_summary)
		colnames(temp_summary) <- gsub( "Pr...t..", "p.value", colnames(temp_summary)  )
		temp_summary$q.value <- p.adjust(temp_summary$p.value, method = "BH")		
		clr_cofound_summary <- rbind(clr_cofound_summary, temp_summary)		
	}	
}
clr_cofound_test_anova <- clr_cofound_test_anova[complete.cases(clr_cofound_test_anova$Pr_Chi),]
clr_cofound_test_anova$q.value <- p.adjust(clr_cofound_test_anova$Pr_Chi,method = "BH")
write.table(file="clr_cofound_test_anova.tsv",clr_cofound_test_anova,sep = "\t", col.names = TRUE, row.names = F)
write.table(file="clr_cofound_summary.tsv",clr_cofound_summary,sep = "\t", col.names = TRUE, row.names = F)

####################################################
#####################   GMPR   #####################

#### Test the individual features using the model and usign the Anova function to check if its significant
gmpr_test_anova <- data.frame()
for(i in  Features){

	subphylo_df <- psmelt( prune_taxa(i,physeq.GMPR) )	
	temp_res <- data.frame()
	for(j in Variables){
		tempList <- nb_glm_function(Y.var = "Abundance", X.var = j, in.df=subphylo_df )
		tempAnova <- data.frame(Feature = i,tempList$anova)
		temp_res <- rbind(temp_res,tempAnova)
		rm(tempList,tempAnova)
	}
	temp_res$q.value <- p.adjust(temp_res$Pr..Chisq., method = "BH")
	gmpr_test_anova <- rbind(gmpr_test_anova,temp_res)
	rm(temp_res,subphylo_df)
}
gmpr_test_anova <- gmpr_test_anova[order(gmpr_test_anova$Pr..Chisq.),]
write.table(file="gmpr_test_anova.tsv",gmpr_test_anova,sep = "\t", col.names = TRUE, row.names = F)
sig_gmpr_test_anova <- gmpr_test_anova[gmpr_test_anova$q.value <= 0.1,]


#### Cofound the analysis ####
gmpr_cofound_test_anova <- data.frame()
gmpr_cofound_summary <- data.frame()
for(i in unique(sig_gmpr_test_anova$Feature)){

	tempsig <- sig_gmpr_test_anova[sig_gmpr_test_anova$Feature == i,]
	subphylo_df <- psmelt( prune_taxa(i,physeq.GMPR) )	
	
	for(j in unique(tempsig$var)){
		tempList <- nb_cofound_function(Y.var = "Abundance", NULL_model = c("Age_at_visit"), 
				Test_model = c("Age_at_visit",j), in.df = subphylo_df)  
		stast_anova <- data.frame(Feature = i,Variable2Test = j, tempList[["Stats"]])
		gmpr_cofound_test_anova <- rbind(gmpr_cofound_test_anova,stast_anova)
	
		temp_summary <- summary(tempList[["Model"]])
		temp_summary <- data.frame(temp_summary$coefficients$cond)
		temp_summary <- data.frame(Feature = i, Variable2Test = j, Varaible = rownames(temp_summary) ,temp_summary)
		colnames(temp_summary) <- gsub( "Pr...z..", "p.value", colnames(temp_summary)  )
		temp_summary$q.value <- p.adjust(temp_summary$p.value, method = "BH")
		
		gmpr_cofound_summary <- rbind(gmpr_cofound_summary, temp_summary)		


	}
	
}
gmpr_cofound_test_anova <- gmpr_cofound_test_anova[complete.cases(gmpr_cofound_test_anova$Pr_Chi),]
gmpr_cofound_test_anova$q.value <- p.adjust(gmpr_cofound_test_anova$Pr_Chi,method = "BH")
write.table(file="gmpr_cofound_test_anova.tsv",gmpr_cofound_test_anova,sep = "\t", col.names = TRUE, row.names = F)
write.table(file="gmpr_cofound_summary.tsv",gmpr_cofound_summary,sep = "\t", col.names = TRUE, row.names = F)


#######################################################################################################################################################
##########################################                   Cofound by the age                         ######################################
#######################################################################################################################################################

#### CLR ####
df2plot <- parser_biomarkers(cofound_summary = clr_cofound_summary, cofound_test_anova = clr_cofound_test_anova, 
			FDR = 0.1, Model_FDR_fil = T, remove_var = c("Age_at_visit","(Intercept)")) 

clr_plot <- ggplot(df2plot, aes(Varaible, Feature, shape = Dominance, fill = Dominance, size = abs(Estimate))) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		geom_point(color = "black") +  # black outline for solid shapes
		scale_shape_manual(values = c("Decrease" = 25, "Increase" = 24)) + 
		scale_fill_manual(values = c("Decrease" = "red", "Increase" = "blue")) + 
		ggtitle("Taxonomic associations")

clr_plot
ggsave(file = "clr_plot.pdf", clr_plot , width = 10, height =  4 )
save(file = "clr_plot.RData", clr_plot )

#### GMPR ####
df2plot <- parser_biomarkers(cofound_summary = gmpr_cofound_summary, cofound_test_anova = gmpr_cofound_test_anova, 
			FDR = 0.1, Model_FDR_fil = T, remove_var = c("Age_at_visit","(Intercept)")) 
gmpr_plot <- ggplot(df2plot, aes(Varaible, Feature, shape = Dominance, fill = Dominance, size = abs(Estimate))) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		geom_point(color = "black") +  # black outline for solid shapes
		scale_shape_manual(values = c("Decrease" = 25, "Increase" = 24)) + 
		scale_fill_manual(values = c("Decrease" = "red", "Increase" = "blue")) + 
		ggtitle("Taxonomic associations")

ggsave(file = "gmpr_plot.pdf", gmpr_plot , width = 10, height =  4 )
save(file = "gmpr_plot.RData", gmpr_plot )

