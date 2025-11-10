set.seed(12345)
library("phyloseq")
library("ggplot2")
library("microbiome")

working_dir <- "~/Postdoc_Raes/Projects/giant_cohort_spa" 

source(paste0(working_dir,"/functions/Extended_Statistical_Toolkit_functions.R"))
source(paste0(working_dir,"/functions/supplementary_figures_functions.R"))

##############################################################################
############################        SCRIPT        ############################
##############################################################################

#load(paste0(working_dir,"/1_infiles/Ileum_biopsies/physeq.ileum.genus.rar.RData"))
#in_phylo <- physeq.ileum.genus.rar
in_phylo <- readRDS(paste0(working_dir,"/1_infiles/Ileum_biopsies/GTDB_r220/physeq_genus.rar.rds"))
Metadata <- sample_data(in_phylo)

### Estimate the alpha diversity ###
Evenness <- evenness(in_phylo, index = "pielou" )
rar_diversity <- estimate_richness(in_phylo, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
rar_diversity$pielou <- Evenness[match(rownames(rar_diversity),rownames(Evenness)),]
Diversity <- rar_diversity[match(rownames(Metadata),rownames(rar_diversity)),]
DiverMeta <- cbind(Metadata,Diversity)
colnames(DiverMeta) <- gsub("Observed","Richness",colnames(DiverMeta))
colnames(DiverMeta) <- gsub("pielou","Pielou",colnames(DiverMeta))

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



