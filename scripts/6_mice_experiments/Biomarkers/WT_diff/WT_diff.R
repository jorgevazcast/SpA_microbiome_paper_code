set.seed(12345)
library(phyloseq)
library(microbiome)

path_project <-  "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

source(paste0(path_project,"/functions/Utils.R"))
source(paste0(path_project,"/functions/beta_diver_functions.R"))
source(paste0(path_project,"/functions/plot_beta_diver_functions.R"))
source(paste0(path_project,"/functions/Extended_Statistical_Toolkit_functions.R"))
source(paste0(path_project,"/functions/supplementary_figures_functions.R"))

###########################################################################################################################
##########################################             FUNCTION                     ######################################
###########################################################################################################################

WTdiff_func <- function(phylo_in,week,Prevalence_filter){
	Metadata <- sample_data(phylo_in)
	Metadata <- Metadata[Metadata$Time == week,] 
	Metadata <- subset(Metadata, enterotype != "not_colonized"  )
	Metadata$is_control <- NULL
	Metadata <- Metadata[complete.cases(Metadata$enterotype),]
	
	#Metadata$Sample_ID <- sapply(Metadata$Sample,function(x){unlist(strsplit(x,"_"))[1]})
	# subdf <- subset(Metadata, Time == "w12" | Time == "w3" )
	B1_prev <- prevalence(prune_samples( rownames(Metadata[Metadata$enterotype == "B1",]) , phylo_in ))
	B2_prev <- prevalence(prune_samples( rownames(Metadata[Metadata$enterotype == "B2",]) , phylo_in ))
	
	B1_prev <- names(B1_prev[B1_prev >= 0.2])
	B2_prev <- names(B2_prev[B2_prev >= 0.2])
		
	in_phylo <- prune_taxa( sort(unique(c(B1_prev,B2_prev))), phylo_in )
	in_phylo <- prune_samples( rownames(Metadata), phylo_in )
	
	phyloDF <- psmelt(in_phylo)

	return_df <- data.frame() 
	for( i in taxa_names(in_phylo)){
		tempdf <- subset( phyloDF, OTU == i )
		res_df <- cont_res_wt(subDF = tempdf,discrete = "enterotype",continuous = "Abundance", Paired = F)
		res_df <- data.frame( OTU = i, res_df  )
		return_df <- rbind(return_df , res_df )
		rm(res_df, tempdf)
	}

	return_df <- return_df[order(return_df$p.value),]
	return_df$q.value <- p.adjust(return_df$p.value,method = "BH")
	return(return_df)
}
		


###########################################################################################################################
##########################################             READ DATA                     ######################################
###########################################################################################################################

EE <- "ee2"
###############################
####  Read the infiles   ######
phyloseq_in <- readRDS(file=paste0(path_project,"/1_infiles/mice_PacBio_data/phyloseq_20k_",EE,"/physeq_sv_decontam_Qfemto.rds"))
taxa_names(phyloseq_in) <- paste0( c(tax_table(phyloseq_in)[,7]) , "_sv_" , 1:ntaxa(phyloseq_in)) 
taxa_names(phyloseq_in) <- gsub(" ","_",taxa_names(phyloseq_in))

Sample2Use <- sample_names(phyloseq_in)
Sample2Use <- Sample2Use[!Sample2Use %in% "Kinnex16S_Fwd_01..Kinnex16S_Rev_13"]
Sample2Use <- Sample2Use[!Sample2Use %in% "X44_batch1"]
phyloseq_in <- prune_samples(Sample2Use,phyloseq_in)

################################################
### Subset only taxa prevalent in the Donors ###
#phylo_donors <- subset_samples(phyloseq_in , Time == "donor" )
#Dprev <- prevalence(phylo_donors)
#Dprev[grep("fragilis",names(Dprev))]
#phyloseq_in <- prune_taxa( names(Dprev[Dprev >= 0.2]) , phyloseq_in )

WT_w12 <- WTdiff_func(phylo_in = phyloseq_in,week = "w12",Prevalence_filter = 0.2)
WT_w3 <- WTdiff_func(phylo_in = phyloseq_in,week = "w3",Prevalence_filter = 0.2)


save(file = "WT_w12.RData" , WT_w12)
save(file = "WT_w3.RData" , WT_w3)

write.table(file= "WT_w12.tsv" ,WT_w12,row.names = F,sep="\t")
write.table(file= "WT_w3.tsv" ,WT_w3,row.names = F,sep="\t")






