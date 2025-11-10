set.seed(12345)
library("phyloseq")
library("ggplot2")
library("microbiome")
library("vegan")
working_dir <- "~/Postdoc_Raes/Projects/giant_cohort_spa" 

source(paste0(working_dir,"/functions/Extended_Statistical_Toolkit_functions.R"))
source(paste0(working_dir,"/functions/supplementary_figures_functions.R"))

load_RData <- function(file) {
      # Load the data into a temporary environment
      temp_env <- new.env()
      load(file, envir = temp_env)
  
      # Get the first object (assuming there's only one)
      obj <- get(ls(temp_env)[1], envir = temp_env)
  
      # Return the renamed object
      return(obj)
}

##############################################################################
############################        SCRIPT        ############################
##############################################################################

# /home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/2_alpha_diversity/Colon/Diversity_test.R
# load(paste0(working_dir,"/1_infiles/Colon_biopsies/physeq.colon.genus.rar.RData"))
in_phylo <- load_RData(paste0(working_dir,"/1_infiles/QMP_GMM/physeq.qmp.gmm.RData"))
# in_phylo <- physeq.colon.genus.rar
Metadata <- sample_data(in_phylo)

# Extract the count/abundance matrix
otu_mat <- as(otu_table(in_phylo), "matrix")

# Transpose if taxa are rows
if (taxa_are_rows(in_phylo)) {
	otu_mat <- t(otu_mat)
}

# Compute alpha diversity indices using vegan
shannon 		<- diversity(otu_mat, index = "shannon")
simpson 		<- diversity(otu_mat, index = "simpson")
inv_simpson 	<- diversity(otu_mat, index = "invsimpson")
observed 		<- specnumber(otu_mat)
#chao1 			<- estimateR(otu_mat)["Chao1", ]

# Estimate evenness using microbiome package
Evenness <- evenness(in_phylo, index = "pielou")

# Merge diversity metrics into a single data frame
rar_diversity <- data.frame(
	Richness		= observed,
	Shannon			= shannon,
	Simpson			= simpson,
	InvSimpson		= inv_simpson,
#	Chao1			= chao1,
	Pielou			= Evenness[match(rownames(otu_mat), rownames(Evenness)), ]
)

# Match metadata with diversity results
Diversity <- rar_diversity[match(rownames(Metadata), rownames(rar_diversity)), ]
DiverMeta <- cbind(Metadata, Diversity)

Diversity_index <- c("Richness","InvSimpson","Pielou")

############################
###### Test Diagnosis ######
res_Diagnosis <- data.frame()
for(i in Diversity_index){
	tempRes <- cont_res_kw(subDF = DiverMeta,discrete = "Diagnosis",continuous = i)
	res_Diagnosis <- rbind( res_Diagnosis , tempRes )

	
}
res_Diagnosis$p.adjust <-  p.adjust(res_Diagnosis$p.value,method="BH")
write.table(res_Diagnosis, file = "KW_Diagnosis_res.tsv", sep = "\t", row.names = F)

##########################
###### Test Disease ######
res_Disease <- data.frame()
for(i in Diversity_index){
	tempRes <- cont_res_wt(subDF = DiverMeta,discrete = "Disease",continuous = i)
	res_Disease <- rbind( res_Disease , tempRes )
}
res_Disease$p.adjust <-  p.adjust(res_Disease$p.value,method="BH")
write.table(res_Disease, file = "WT_Disease_res.tsv", sep = "\t", row.names = F)

##########################
###### Test Disease_activity ######
res_Disease_activity <- data.frame()
for(i in Diversity_index){
	tempRes <- cont_res_kw(subDF = DiverMeta,discrete = "Disease_activity",continuous = i)
	res_Disease_activity <- rbind( res_Disease_activity , tempRes )

	
}
res_Disease_activity$p.adjust <-  p.adjust(res_Disease_activity$p.value,method="BH")
write.table(res_Disease_activity, file = "KW_Disease_activity.tsv", sep = "\t", row.names = F)



####################################################
############### Plot and poshoc test ###############

#######################
###### Diagnosis ######
colors2use <- color_varaibles("Diagnosis")
OrderFeat = names(colors2use)
ret_statsplot <- plot_stats_kw(DiverMeta = DiverMeta, colors2use=colors2use, stats_df=res_Diagnosis, Var= "Diagnosis", OrderFeat = names(colors2use) )

write.table(ret_statsplot$dunnTest_res, file = "dunnTest_Diagnosis_res.tsv", sep = "\t", row.names = F)

list_plots <- ret_statsplot[["list_plots"]]
save(file = "plots_Diagnosis.RData",list_plots)

plots_Diagnosis <-ggarrange(plotlist = list_plots, nrow = 1, common.legend = T, legend = "bottom")
ggsave(file="plots_Diagnosis.pdf",plots_Diagnosis,width = 12, height = 5 )

##############################
###### Disease_activity ######
colors2use <- color_varaibles("Disease_activity")
OrderFeat = names(colors2use)
ret_statsplot <- plot_stats_kw(DiverMeta = DiverMeta, colors2use=colors2use, stats_df=res_Disease_activity, Var= "Disease_activity", OrderFeat = names(colors2use) )

write.table(ret_statsplot$dunnTest_res, file = "dunnTest_Disease_activity_res.tsv", sep = "\t", row.names = F)

list_plots <- ret_statsplot[["list_plots"]]
save(file = "plots_Disease_activity.RData",list_plots)

plots_Disease_activity <-ggarrange(plotlist = list_plots, nrow = 1, common.legend = T, legend = "bottom")
ggsave(file="plots_Disease_activity.pdf",plots_Disease_activity,width = 12, height = 5 )

#####################
###### Disease ######
colors2use <- color_varaibles("Disease")
OrderFeat = names(colors2use)
list_plots <- plot_stats_wt(DiverMeta = DiverMeta, colors2use=colors2use, stats_df=res_Disease, Var= "Disease", OrderFeat = names(colors2use) )

save(file = "plots_Disease.RData",list_plots)
plots_Disease <-ggarrange(plotlist = list_plots, nrow = 1, common.legend = T, legend = "bottom")
ggsave(file="plots_Disease.pdf",plots_Disease,width = 12, height = 5 )



