#### Perform the DA analysis and return the significant taxa #####

set.seed(12345)
library(microbiome)
library(ACAT) # ACAT takes a series of p-values ​​(e.g., from different statistical tests applied to the same hypothesis or gene) and combines them into a single statistic using the Cauchy distribution.

argscomd = commandArgs(trailingOnly=TRUE)
print(argscomd[1])

GMPR.file <- as.character(argscomd[1])
rar.file <- as.character(argscomd[2])
above10000.file <- as.character(argscomd[3])
clr.file <- as.character(argscomd[4])
Variable = as.character(argscomd[5]) ## Variable = "Disease_activity" # Histological_Inflammation Histology Disease
Tax = as.character(argscomd[6]) ## 

cat("\n Parameters: ", "\n", GMPR.file, "\n", rar.file, "\n", above10000.file, "\n", clr.file,"\n", Variable,"\n", Tax,"\n" )

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


### Create the data.frames
List_Tables <- filter_and_match_data_phyloseq(
		physeq.GMPR = physeq.GMPR,
		physeq.rar = physeq.rar,
		physeq.all = physeq.all,
		physeq.clr = physeq.clr, 
		inMetadata = data.frame(sample_data(physeq.rar)), 
		Var2Use = Variable, 
		prev_cutoff =  0.2
	)

taxa_GMPR <- List_Tables[["taxa_GMPR"]]
taxa_raw <- List_Tables[["taxa_raw"]]
taxa_rar <- List_Tables[["taxa_rar"]]
taxa_clr <- List_Tables[["taxa_clr"]] ### The CLR might have less taxa
Metadata <- List_Tables[["Metadata"]]

###################################
########    Filter tests   ########

Ncats <- length(unique(Metadata[[Variable]]))
if(Ncats == 2){
	Test <- c("ALDEx2","DESeq2","GMPR_wt","RAR_wt","CLR_ttest","CLR_wt","ANCOMBC")
}else if(Ncats > 2){
	Test <- c("ALDEx2","DESeq2","GMPR_kw","RAR_kw","CLR_anova","CLR_kw","ANCOMBC")
}
remove_test <- "ALDEx2"

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
		Test = c("ALDEx2","DESeq2","GMPR_wt","RAR_wt","CLR_ttest","CLR_wt","ANCOMBC"),
		Meta = Metadata,
		in_phylo = physeq.rar,
		taxa.rar = taxa_rar,
		taxa.GMPR = taxa_GMPR,
		taxa.clr = taxa_clr,
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

Plot_summary_GMPR <- plot_significant_features( in.meta = Metadata, in.table = taxa_GMPR, pval_table = q_values_table_sig, 
	Var2Use = Variable, cutoff = p_val_cut_off, colors=colors_boxplot , log_data = T )
save(file = "Plot_summary_GMPR.RData", Plot_summary_GMPR )
ggsave(file = "Plot_summary_GMPR.pdf", Plot_summary_GMPR, width = 16, height =4)

Plot_summary_CLR <- plot_significant_features( in.meta = Metadata, in.table = taxa_clr, pval_table = q_values_table_sig, 
	Var2Use = Variable, cutoff = p_val_cut_off,colors=colors_boxplot  )
save(file = "Plot_summary_CLR.RData", Plot_summary_CLR )
ggsave(file = "Plot_summary_CLR.pdf", Plot_summary_CLR, width = 16, height = 4)




