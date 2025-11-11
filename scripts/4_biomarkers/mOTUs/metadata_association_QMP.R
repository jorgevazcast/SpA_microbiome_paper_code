set.seed(12345)
library(glmmTMB)    
library(glmnet)
library(lmerTest)
library(microbiome)

argscomd = commandArgs(trailingOnly=TRUE)

rar.file <- as.character(argscomd[1])
Variable = as.character(argscomd[2]) ## Variable = "Disease_activity" # Histological_Inflammation Histology Variable = "Disease"
sig_features =  as.character(argscomd[3]) 

cat("\n Parameters: ", "\n", rar.file,"\n", Variable,"\n", sig_features,"\n" )

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
physeq.rar <- read_input_phyloseq(rar.file)

taxa_names(physeq.rar) <- gsub("[^[:alnum:]]", "_", taxa_names(physeq.rar))
taxa_names(physeq.rar) <- gsub(" ", "_", taxa_names(physeq.rar))	

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
load("~/github_shared_code_and_publications/SpA_microbiome_paper_code/Metadata/metagenomic_varaibles_associations/Shotgun_varaibles/Variables2Use.RData")
load("~/github_shared_code_and_publications/SpA_microbiome_paper_code/Metadata/metagenomic_varaibles_associations/Shotgun_varaibles/Var2eliminate.RData")
var2use <- as.character(Variables2Use$Var)
Variables2eliminate <- Var2eliminate$Var
Variables <- var2use[!var2use %in% Variables2eliminate]

Features <- unique(q_values_table_sig$Feature)

####################################################
#####################   RAR   #####################

#### Test the individual features using the model and usign the Anova function to check if its significant
rar_test_anova <- data.frame()
for(i in  Features){

	subphylo_df <- psmelt( prune_taxa(i,physeq.rar) )	
	temp_res <- data.frame()
	for(j in Variables){
		tempList <- nb_glm_function(Y.var = "Abundance", X.var = j, in.df=subphylo_df )
		tempAnova <- data.frame(Feature = i,tempList$anova)
		temp_res <- rbind(temp_res,tempAnova)
		rm(tempList,tempAnova)
	}
	temp_res$q.value <- p.adjust(temp_res$Pr..Chisq., method = "BH")
	rar_test_anova <- rbind(rar_test_anova,temp_res)
	rm(temp_res,subphylo_df)
}
rar_test_anova <- rar_test_anova[order(rar_test_anova$Pr..Chisq.),]
write.table(file="rar_test_anova.tsv",rar_test_anova,sep = "\t", col.names = TRUE, row.names = F)

sig_rar_test_anova <- rar_test_anova[rar_test_anova$q.value <= 0.1,]
sig_rar_test_anova <- sig_rar_test_anova[complete.cases(sig_rar_test_anova),]


#### Cofound the analysis ####
rar_cofound_test_anova <- data.frame()
rar_cofound_summary <- data.frame()
for(i in unique(sig_rar_test_anova$Feature)){

	tempsig <- sig_rar_test_anova[sig_rar_test_anova$Feature == i,]
	subphylo_df <- psmelt( prune_taxa(i,physeq.rar) )	
	
	for(j in unique(tempsig$var)){
	
		tempList <- nb_cofound_function(Y.var = "Abundance", NULL_model = c("Water","BMI","Sex"), 
				Test_model = c("Water","BMI","Sex",j), in.df = subphylo_df)  
				 
		stast_anova <- data.frame(Feature = i,Variable2Test = j, tempList[["Stats"]])
		rar_cofound_test_anova <- rbind(rar_cofound_test_anova,stast_anova)
	
		temp_summary <- summary(tempList[["Model"]])
		temp_summary <- data.frame(temp_summary$coefficients$cond)
		temp_summary <- data.frame(Feature = i, Variable2Test = j, Varaible = rownames(temp_summary) ,temp_summary)
		colnames(temp_summary) <- gsub( "Pr...z..", "p.value", colnames(temp_summary)  )
		temp_summary$q.value <- p.adjust(temp_summary$p.value, method = "BH")
		
		rar_cofound_summary <- rbind(rar_cofound_summary, temp_summary)		
		rm(temp_summary,tempList,stast_anova)

	}
	
}
rar_cofound_test_anova <- rar_cofound_test_anova[complete.cases(rar_cofound_test_anova$Pr_Chi),]
rar_cofound_test_anova$q.value <- p.adjust(rar_cofound_test_anova$Pr_Chi,method = "BH")
write.table(file="rar_cofound_test_anova.tsv",rar_cofound_test_anova,sep = "\t", col.names = TRUE, row.names = F)
write.table(file="rar_cofound_summary.tsv",rar_cofound_summary,sep = "\t", col.names = TRUE, row.names = F)


#######################################################################################################################################################
##########################################                   Cofound by the age                         ######################################
#######################################################################################################################################################

#### RAR ####
df2plot <- parser_biomarkers(cofound_summary = rar_cofound_summary, cofound_test_anova = rar_cofound_test_anova, 
			FDR = 0.1, Model_FDR_fil = T, remove_var = c("Water","BMI","SexMale","(Intercept)")) 
rar_plot <- ggplot(df2plot, aes(Varaible, Feature, shape = Dominance, fill = Dominance, size = abs(Estimate))) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		geom_point(color = "black") +  # black outline for solid shapes
		scale_shape_manual(values = c("Decrease" = 25, "Increase" = 24)) + 
		scale_fill_manual(values = c("Decrease" = "red", "Increase" = "blue")) + 
		ggtitle("Taxonomic associations")
rar_plot
ggsave(file = "rar_plot.pdf", rar_plot , width = 10, height =  6 )
save(file = "rar_plot.RData", rar_plot )
save(file = "df2plot.RData", df2plot )



