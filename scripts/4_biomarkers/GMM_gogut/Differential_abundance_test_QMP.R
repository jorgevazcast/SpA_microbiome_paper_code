
#### Perform the DA analysis and return the significant taxa #####

set.seed(12345)
library(microbiome)
library(ACAT) # ACAT takes a series of p-values ​​(e.g., from different statistical tests applied to the same hypothesis or gene) and combines them into a single statistic using the Cauchy distribution.

argscomd = commandArgs(trailingOnly=TRUE)

rar.file <- as.character(argscomd[1]) 
Variable = as.character(argscomd[2])
cat("\n Parameters: ", "\n", rar.file,"\n", Variable,"\n", "\n" )

#################################
######## q-value cut-off ########
p_val_cut_off <- 0.1

#######################################################################################################################################
##########################################              Load the functions                       ######################################
#######################################################################################################################################

working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

source(paste0(working_dir,"/Functions/data_processing_functions.R"))
source(paste0(working_dir,"/Functions/statistical_test_functions.R"))
source(paste0(working_dir,"/functions/CATDB_supplementary_functions.R"))


#######################################################################################################################################
##########################################                  Read the data                        ######################################
#######################################################################################################################################

### Read the phyloseq objects
physeq.rar <- read_input_phyloseq(rar.file)
MetadataIn <- data.frame(sample_data(physeq.rar))
MetadataIn$Response <- ifelse(MetadataIn$Treatment_response == "Non_responder" , "Non_responder", "Responder" )

fake_taxa <- data.frame( Dummy = 1, Names = taxa_names(physeq.rar) )
rownames(fake_taxa) <-  taxa_names(physeq.rar)

physeq.rar <- create_phyloseq_object_fun( data = as.data.frame( as.matrix(physeq.rar@otu_table) ) , 
				tax = fake_taxa, sample_dat = MetadataIn  )
### Create the data.frames
List_Tables <- filter_and_match_data_phyloseq_QMP(
		physeq.rar = physeq.rar,
		inMetadata = data.frame(sample_data(physeq.rar)), 
		Var2Use = Variable, 
		prev_cutoff =  0.2
	)
taxa_rar <- List_Tables[["taxa_rar"]]
Metadata <- List_Tables[["Metadata"]]

###################################
########    Filter tests   ########

Ncats <- length(unique(Metadata[[Variable]]))
if(Ncats == 2){
	Test <- c("ALDEx2","DESeq2","RAR_wt","ANCOMBC")
}else if(Ncats > 2){
	Test <- c("ALDEx2","DESeq2","RAR_kw","ANCOMBC")
}
remove_test <- c("ALDEx2","DESeq2")
Test <- Test[!Test %in% remove_test]
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
		in_phylo = physeq.rar,
		taxa.rar = taxa_rar,
#		taxa.GMPR = taxa_GMPR,
#		taxa.clr = taxa_clr,
		Var = Variable,
		QMPdata = T
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
q_values_table <- q_values_table[complete.cases(q_values_table),]
q_values_table <- subset(q_values_table, ACAT <= p_val_cut_off ) ### The ACAT must be significant
Sig_taxas <- apply(   q_values_table[, 2:(ncol(q_values_table)-1)], 1, function(x) { any(x <= p_val_cut_off)  } ) ### Other test must be significat
q_values_table_sig <- q_values_table[Sig_taxas,]

Plot_summary_RAR <- plot_significant_features( in.meta = Metadata, in.table = taxa_rar, pval_table = q_values_table_sig, 
	Var2Use = Variable, cutoff = p_val_cut_off, colors=colors_boxplot , log_data = T )
save(file = "Plot_summary_RAR.RData", Plot_summary_RAR )
ggsave(file = "Plot_summary_RAR.pdf", Plot_summary_RAR, width = 16, height = 4)


