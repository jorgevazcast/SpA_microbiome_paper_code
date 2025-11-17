set.seed(12345)
library(phyloseq)
library(microbiome)


qval_cutoff  <- 0.1
ESize_r_cutoff <- 0.3
EE <- "ee2"

path_project <-  "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

###############################
####  Read the infiles   ######

source(paste0(path_project,"/functions/Extended_Statistical_Toolkit_functions.R"))
source(paste0(path_project,"/functions/supplementary_figures_functions.R"))

phyloseq_in <- readRDS(file=paste0(path_project,"/1_infiles/mice_PacBio_data/phyloseq_20k_",EE,"/physeq_sv_decontam_Qfemto.rds"))
taxa_names(phyloseq_in) <- paste0( c(tax_table(phyloseq_in)[,7]) , "_sv_" , 1:ntaxa(phyloseq_in)) 
taxa_names(phyloseq_in) <- gsub(" ","_",taxa_names(phyloseq_in))

Sample2Use <- sample_names(phyloseq_in)
Sample2Use <- Sample2Use[!Sample2Use %in% "Kinnex16S_Fwd_01..Kinnex16S_Rev_13"]
Sample2Use <- Sample2Use[!Sample2Use %in% "X44_batch1"]
phyloseq_in <- prune_samples(Sample2Use,phyloseq_in)

load(paste0(path_project,"/6_mice_experiments/Biomarkers/WT_diff/WT_w12.RData"))

WT_w12 <- WT_w12[complete.cases(WT_w12$q.value),]
sig_features <- WT_w12[WT_w12$q.value < qval_cutoff,]
sig_features <- sig_features[sig_features$EffectSize_r > ESize_r_cutoff,]

sub_phylo <- prune_taxa( sig_features$OTU , phyloseq_in)
sub_phylo <- subset_samples( sub_phylo , Time == "w12"  )


Exclude_variables <- c("Qfemto")
Metadata <- data.frame(sample_data( sub_phylo ))

#### Metadata #####
N_values_NA_cat <- N_categories_NA_values_function(in.meta = Metadata)

### High NA data
#High_NA_var <- N_values_NA_cat[N_values_NA_cat$NA_percentage >= 20,]$Variable

### Single categorie variable
Single_categorie <- subset(N_values_NA_cat,  N_values == 1 )$Variable

### Remove binary variables with less than 1% of imbalance data
min_num_samp_per_cat <- round(nrow(Metadata) * 0.01)
Non_Info_var <- subset(N_values_NA_cat, Var_type == "character" & N_min_categories < min_num_samp_per_cat & N_values == 2 )$Variable

#### Exclude the varaibles ####
Exclude_variables <- unique( c(  Exclude_variables ,  Single_categorie, Non_Info_var ) )
Metadata <- data.frame(Metadata[,!colnames(Metadata) %in% Exclude_variables])

#### Cont variales 
class_cat <- N_categories_NA_values_function(in.meta = Metadata)
Var2Corr <- class_cat[class_cat$Var_type == "numeric",]$Variable

###################################################
#### Correlate sig taxa with all the varaibles ####
df_phylo <- psmelt(sub_phylo)
cor_taxa_metadata <- data.frame()
for(tax in unique(df_phylo$OTU)){
	# tax <- "s_Bacteroides_fragilis_sv_150"
	subdf <- subset( df_phylo, OTU == tax )
	temp_res <- data.frame()
	for(var in Var2Corr){
		# var <-  "synov_CCL2"
		tempcor <- cor_func(cont_var1 = "Abundance",cont_var2 = var, Method = "spearman", subDF = subdf)
		temp_res <- rbind(temp_res, data.frame( OTU = tax, tempcor))
		rm(tempcor)
	}
	temp_res$q.val <- p.adjust(temp_res$p.val, method="BH")
	cor_taxa_metadata <- rbind(cor_taxa_metadata , temp_res)
	rm(tempcor,temp_res,subdf)	
	
}

########################
#### Save the files ####
cor_taxa_metadata <- cor_taxa_metadata[order(cor_taxa_metadata$q.val),]
cor_taxa_metadata$Sig <- ""
cor_taxa_metadata$q.val[is.na(cor_taxa_metadata$q.val)] <- 1  # Poner 1 a los NA
cor_taxa_metadata$Sig <- ifelse(cor_taxa_metadata$q.val < 0.01, "***",
                          ifelse(cor_taxa_metadata$q.val < 0.05, "**",
                          ifelse(cor_taxa_metadata$q.val < 0.1, "*", "")))


save(file = "cor_taxa_metadata.RData" , cor_taxa_metadata)
write.table(file= "cor_taxa_metadata.tsv" ,cor_taxa_metadata,row.names = F,sep="\t")

########################
#### Plot the results ####
cor_taxa_metadata_sig <- subset(cor_taxa_metadata , q.val < qval_cutoff)
cor_taxa_metadata_sig <- cor_taxa_metadata_sig[order(cor_taxa_metadata_sig$estimate ,decreasing = T),]

Cols <- names(sort(table(cor_taxa_metadata_sig$Var2)))
Rows <- names(sort(table(cor_taxa_metadata_sig$OTU)))

RhoMatrix <- matrix( NA , length(Rows) , length(Cols)  )
colnames(RhoMatrix) <- Cols
rownames(RhoMatrix) <- Rows



sig_features$EffectSize_r <- ifelse( sig_features$Dominant == "B2", (sig_features$EffectSize_r * -1), sig_features$EffectSize_r )

sig_features <- sig_features[order(sig_features$EffectSize_r,decreasing = T),]
sig_features <- sig_features[sig_features$OTU %in% rownames(RhoMatrix), ]
EffectSize_r <- sig_features$EffectSize_r
names(EffectSize_r) <- sig_features$OTU

RhoMatrix <- RhoMatrix[match( names(EffectSize_r) ,  rownames(RhoMatrix) ),]
QvalMatrix <- RhoMatrix

for(i in rownames( RhoMatrix )){
	for(j in colnames( RhoMatrix )){
		RhoMatrix[i,j] <- subset( cor_taxa_metadata, OTU == i & Var2 == j )$estimate
	
		QvalMatrix[i,j] <- subset( cor_taxa_metadata, OTU == i & Var2 == j )$Sig
		
	}

}

annotation_row = data.frame( Dominant = sig_features$Dominant )
ann_colors = list(  Dominant = c(B2 = "#ba1015", B1 = "#f28118" ) )


library(ComplexHeatmap)

p1 <- pheatmap( RhoMatrix,
	#display_numbers = QvalMatrix,
	cluster_rows = FALSE,
	cluster_cols = FALSE,
	row_names_side = "left",
	color = colorRampPalette( rev( c("navy", "white", "firebrick3")) )(50),
	fontsize_number = 17,
	heatmap_legend_param = list( 	title = "Spearman Rho",  at = c(-1, 0, 1), color_bar = "continuous"),
	cell_fun = function(j, i, x, y, width, height, fill) {
	  grid.text(QvalMatrix[i, j], x, y, gp = gpar(fontsize = 12))
	},
	annotation_row = annotation_row,
	annotation_colors = ann_colors

)


p2 <- rowAnnotation(EffectSize_r = anno_barplot( EffectSize_r  ) )

pdf("Heatmap.pdf", width = 13 , height = 22)
	draw(p1 + p2 , heatmap_legend_side = "right", column_title = "Correlation between Significant Taxa and Metadata")

dev.off()


##################################################################################################
##################################################################################################

cor_taxa_metadata_sub <- subset(cor_taxa_metadata, Var2 == "colon_TNF")
cor_taxa_metadata_sub <- subset(cor_taxa_metadata_sub, q.val < qval_cutoff)
cor_taxa_metadata_sub$colon_TNF <- ifelse( cor_taxa_metadata_sub$estimate > 0, "Postive", "Negative" )
cor_taxa_metadata_sub <- cor_taxa_metadata_sub[order(cor_taxa_metadata_sub$estimate),]
cor_taxa_metadata_sub$OTU <- factor( as.character(cor_taxa_metadata_sub$OTU) , cor_taxa_metadata_sub$OTU )

pImportance <- ggplot(cor_taxa_metadata_sub, aes(x = OTU, y = estimate, fill = colon_TNF)) + 
	geom_bar(stat = "identity") +
	#geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
	coord_flip() + 
	theme_bw() + 
	ylab("Spearman Rho") +
	ggtitle(paste("w12", "Taxa associated with colon TNF levels")) +
	scale_fill_manual(values = c("red","blue")) #+ theme_classic()

pImportance	
ggsave( file = "Taxa_correlation.pdf" , pImportance, width = 9, height = 17 )	
saveRDS(file=paste0("w12",".cor_taxa_metadata.plot.rds"),pImportance)







