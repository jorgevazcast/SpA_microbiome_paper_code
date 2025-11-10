
#### Perform the DA analysis and return the significant taxa #####

set.seed(12345)
library(microbiome)
library(ACAT) # ACAT takes a series of p-values ​​(e.g., from different statistical tests applied to the same hypothesis or gene) and combines them into a single statistic using the Cauchy distribution.

argscomd = commandArgs(trailingOnly=TRUE)

log.file <- as.character(argscomd[1])
raw.file <- as.character(argscomd[2])
Variable = as.character(argscomd[3]) ##  Variable = "Disease"  # Variable = "Disease_activity" # Histological_Inflammation Histology
#log.file <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/Metabolomics/physeq.metabolites.norm.RData"
#raw.file <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/Metabolomics/physeq.metabolites.RData"


#cat("\n Parameters: ", "\n", GMPR.file, "\n", rar.file, "\n", above10000.file, "\n", clr.file,"\n", Variable,"\n", Tax,"\n" )
cat("\n Parameters: ", "\n",log.file, "\n", raw.file,"\n", Variable,"\n", "\n" )
#################################
######## q-value cut-off ########
p_val_cut_off <- 0.1

#######################################################################################################################################
##########################################              Load the functions                       ######################################
#######################################################################################################################################

path_func <- "/home/luna.kuleuven.be/u0141268/github_projects/CATBD"
#source(paste0(path_func,"/Functions/functions_Biomarker.R"))
source(paste0(path_func,"/Functions/data_processing_functions.R"))
source(paste0(path_func,"/Functions/statistical_test_functions.R"))

path_func_SpA <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa"
source(paste0(path_func_SpA,"/functions/CATDB_supplementary_functions.R"))


#######################################################################################################################################
##########################################                  Read the data                        ######################################
#######################################################################################################################################

### Read the phyloseq objects
physeq.log <- read_input_phyloseq(log.file)
physeq.raw <- read_input_phyloseq(raw.file)
MetadataIn <- data.frame(sample_data(physeq.log))

### Create the data.frames
List_Tables <- filter_and_match_data_phyloseq_metabolome(
		physeq.raw = physeq.raw,
		physeq.norm = physeq.log,		
		inMetadata = data.frame(sample_data(physeq.raw)), 
		Var2Use = Variable, 
		prev_cutoff =  0.2
	)
taxa_norm <- List_Tables[["taxa_norm"]]
taxa_raw <- List_Tables[["taxa_raw"]]
Metadata <- List_Tables[["Metadata"]]

###################################
########    Filter tests   ########

Ncats <- length(unique(Metadata[[Variable]]))
if(Ncats == 2){
	Test <- c("ANCOMBC" , "log_ttest", "log_wt", "raw_wt" )
}else if(Ncats > 2){
	Test <- c("ANCOMBC" , "log_anova", "log_kw", "raw_kw" )
}
remove_test <- c("raw_wt", "ANCOMBC" )

#########################################################
##########    Putative contamination taxa    ############
#########################################################

contamination <- contamination_taxa()

#########################################################
########## Setup the pairwise comparassions  ############
#########################################################

list_colors <- Disease_color_function(metadata_variable = Variable,Metadata = Metadata)
colors_boxplot <- list_colors[["colors_boxplot"]]
levels_boxplot <- list_colors[["levels_boxplot"]]

#######################################################################################################################################
##########################################            Create the outfiles directory              ######################################
#######################################################################################################################################

dir2create <- paste0(getwd(),"/",Variable)
dir.create(file.path(dir2create ), showWarnings = FALSE)
setwd(file.path(dir2create))

####################################################################################################################################################
##########################################                  Statistical abundance tests                        #####################################
####################################################################################################################################################

Test_res <- RUN_tests(  
		Test = Test,
		Meta = Metadata,
		in_phylo = physeq.raw,
		taxa.log = taxa_norm,
		taxa.raw = taxa_raw,
		Var = Variable
)

####################################################################################################################################################
##########################################                     Join p-value table                              #####################################
####################################################################################################################################################

if(!is.null(remove_test)){  Test_res <-  Test_res[!names(Test_res) %in% remove_test]   }

pvalues_list <- join_pvalues_function(Test_list = Test_res, Ncat =  Ncats )

p_values_table <- pvalues_list$p_values_table_res

q_values_table <- pvalues_list$q_values_table_res

####################################################################################################################################################
##########################################                 Plot the significat results                         #####################################
####################################################################################################################################################

### Subset the signifcat results
q_values_table <- subset(q_values_table, ACAT <= p_val_cut_off ) ### The ACAT must be significant
Sig_taxas <- apply(   q_values_table[, 2:(ncol(q_values_table)-1)], 1, function(x) { any(x <= p_val_cut_off)  } ) ### Other test must be significat
q_values_table_sig <- q_values_table[Sig_taxas,]

Plot_summary_log <- plot_significant_features( in.meta = Metadata, in.table = taxa_norm, pval_table = q_values_table_sig, 
	Var2Use = Variable, cutoff = p_val_cut_off, colors=colors_boxplot , log_data = F )
save(file = "Plot_summary_log.RData", Plot_summary_log )
ggsave(file = "Plot_summary_log.pdf", Plot_summary_log, width = 16, height = 4)


