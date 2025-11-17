set.seed(12345)
library(microbiome)
library(phyloseq)
library(glmnet)
library(lmerTest)


Disease_color_function <- function(metadata_variable,Metadata){
	if(metadata_variable == "Disease_activity"){
		colors_var = c("blue","red", "darkgreen")
		names(colors_var) <- c("HC" , "SpAHigh","SpALow")
		colors_boxplot=colors_var
		levels_boxplot= unique(sort(as.character(Metadata$Disease_activity)))
		print(metadata_variable)
		print(levels_boxplot)
		
	}

	if(metadata_variable == "Disease"){
		colors_var = c("blue","red")
		names(colors_var) <- c("HC" , "SpA")
		colors_boxplot=colors_var
		levels_boxplot= unique(sort(as.character(Metadata$Disease)))
		print(metadata_variable)
		print(levels_boxplot)
	}

	if(metadata_variable == "Response"){
		colors_var = c("blue","red")
		names(colors_var) <- c("Responder" , "Non_responder")
		colors_boxplot=colors_var
		levels_boxplot= unique(sort(as.character(Metadata$Response)))
		print(metadata_variable)
		print(levels_boxplot)
	}
	list_ret <- list(colors_boxplot,levels_boxplot)
	names(list_ret) <- c("colors_boxplot","levels_boxplot")
	return(list_ret)
}



ASV_naming_function <- function( Sp_annotations, Taxonomy){ # Taxonomy;: "GreenGenes2" "GTDB_r86" "GTDB_r202" "GTDB_r207" "GTDB_r220" "rdp_16" "rdp_19" "silva_v138_2" "Consensus_sp"

	DF_annotation <- data.frame(ASV = rownames(Sp_annotations) , Annotation = Sp_annotations[[Taxonomy]] )
	DF_annotation$ID <- paste0(DF_annotation$Annotation,"_asv_",1:nrow(DF_annotation))
	DF_annotation$ID <- gsub("[^[:alnum:]]", "_", DF_annotation$ID)
	DF_annotation$ID <- gsub(" ", "_", DF_annotation$ID)	
	DF_annotation[7,]
	return(DF_annotation)
}

phyloseq_asv_annotation <- function(in.phylo , Annotation){

	Annotation <- Annotation[match(taxa_names(in.phylo) ,  Annotation$ASV),]
	taxa_names(in.phylo) <- Annotation[["ID"]]
	taxa_names(in.phylo) <- gsub("[^[:alnum:]]", "_", taxa_names(in.phylo))
	taxa_names(in.phylo) <- gsub(" ", "_", taxa_names(in.phylo))	
	return(in.phylo)
}

parser_biomarkers <- function(cofound_summary = data.frame(), cofound_test_anova = data.frame(), FDR = 0.1, Model_FDR_fil = T, remove_var = c("(Intercept)")) {
	
	# Create a unique identifier for each row in the ANOVA test data
	cofound_test_anova$ID <- paste(cofound_test_anova$Feature, cofound_test_anova$Variable2Test)

	# Filter the ANOVA test results by FDR threshold
	sig_cofound_test_anova <- cofound_test_anova[cofound_test_anova$q.value <= FDR, ]

	# Recompute ID column for filtered results (optional redundancy)
	sig_cofound_test_anova$ID <- paste(sig_cofound_test_anova$Feature, sig_cofound_test_anova$Variable2Test)

	# Create ID column for the summary table
	cofound_summary$ID <- paste(cofound_summary$Feature, cofound_summary$Variable2Test)

	# Select only rows in the summary table that passed the significance threshold
	sig_cofound_summary <- cofound_summary[cofound_summary$ID %in% sig_cofound_test_anova$ID, ]

	# Reset row names
	rownames(sig_cofound_summary) <- NULL

	# Remove intercept terms from results
	sig_cofound_summary <- sig_cofound_summary[!sig_cofound_summary$Varaible %in% remove_var, ]

	if(Model_FDR_fil == T){ sig_cofound_summary <- sig_cofound_summary[sig_cofound_summary$q.value <= FDR, ]  }

	# Add a Dominance column to indicate direction of the effect
	sig_cofound_summary$Dominance <- ifelse(sig_cofound_summary$Estimate > 0, "Increase", "Decrease")
	# Sort features by frequency (least to most frequent)
	SortFeatures <- names(sort(table(sig_cofound_summary$Feature), decreasing = FALSE))
	# Sort variables by frequency (most to least frequent)
	SortCat <- names(sort(table(sig_cofound_summary$Varaible), decreasing = TRUE))
	# Prioritize "DiseaseSpA" to appear first in the variable levels
	SortCat <- c("DiseaseSpA", SortCat[!grepl("DiseaseSpA", SortCat)])
	# Set factor levels for plotting/ordering purposes
	sig_cofound_summary$Varaible <- factor(as.character(sig_cofound_summary$Varaible), levels = SortCat)
	sig_cofound_summary$Feature <- factor(as.character(sig_cofound_summary$Feature), levels = SortFeatures)
	# Return the filtered and formatted summary
	return(sig_cofound_summary)
}



