set.seed(12345)
library(phyloseq)
library(microbiome)
library(ComplexHeatmap)


qval_cutoff  <- 0.1
ESize_r_cutoff <- 0.3

load(file = "cor_taxa_metadata.RData")
sigfiles <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/6_mice_experiments/Biomarkers"
load(paste0(sigfiles,"/all_samples/WT_diff/WT_w12.RData"))

WT_w12 <- WT_w12[complete.cases(WT_w12$q.value),]
sig_features <- WT_w12[WT_w12$q.value < qval_cutoff,]
sig_features <- sig_features[sig_features$EffectSize_r > ESize_r_cutoff,]


########################
#### Plot the results ####
cor_taxa_metadata_sig <- subset(cor_taxa_metadata , q.val < qval_cutoff)
cor_taxa_metadata_sig <- cor_taxa_metadata_sig[order(cor_taxa_metadata_sig$estimate ,decreasing = T),]
sigMoreThanVar <- table(cor_taxa_metadata_sig$OTU)
sigMoreThanVar <- names(sigMoreThanVar[sigMoreThanVar > 1])
sigNonTNF <- cor_taxa_metadata_sig[cor_taxa_metadata_sig$Var2 != "colon_TNF",]$OTU

SigVar <- unique(c(sigMoreThanVar , sigNonTNF))

cor_taxa_metadata_sig <- cor_taxa_metadata_sig[cor_taxa_metadata_sig$OTU %in% SigVar,]

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
pdf("Heatmap.subset_var.pdf", width = 13 , height = 9)
	draw(p1 + p2 , heatmap_legend_side = "right", column_title = "Correlation between Significant Taxa and Metadata")
dev.off()







