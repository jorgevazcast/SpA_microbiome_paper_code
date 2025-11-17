set.seed(12345)
library(ggpubr)
library(ggplot2)
library(scales)


#######################################################################################################################################
##########################################                     FUNCTIONS                          #####################################
#######################################################################################################################################


#Value of r	Effect size interpretation
#0.10	Small effect
#0.30	Medium effect
#0.50	Large effect


WTdiff_filter_plot <- function( week, sig_features){

	sig_features <- sig_features[complete.cases(sig_features$OTU),]
	sig_features <- sig_features[complete.cases(sig_features$q.value),]
	sig_features <- sig_features[complete.cases(sig_features$EffectSize_r),]	
		
	sig_features <- sig_features[sig_features$q.value < 0.1,]
	sig_features <- sig_features[sig_features$EffectSize_r > 0.3,]

	sig_features$Effect_size <- ifelse(sig_features$EffectSize_r > 0.3, "Large effect", "Medium effect")
	sig_features$Effect_size <- ifelse(sig_features$EffectSize_r > 0.5, "Large effect", "Medium effect")	
	sig_features$ColorGroup <- paste(sig_features$Dominant, sig_features$Effect_size)	
	sig_features$EffectSize_r <- ifelse(sig_features$Dominant == "B1", sig_features$EffectSize_r, sig_features$EffectSize_r * -1)	
	
	my_colors <- c(
	  "B1 Large effect" = "#f28118",                       # Orange
	  "B1 Medium effect" = scales::alpha("#f28118", 0.5),          # Same color, 50% transparency
	  "B2 Large effect" = "#ba1015",                       # Red-brown
	  "B2 Medium effect" = scales::alpha("#ba1015", 0.5)           # Same color, 50% transparency
	)

	sig_features <- sig_features[order(sig_features$EffectSize_r),]
	sig_features$OTU <- factor( as.character(sig_features$OTU) , sig_features$OTU )
	pImportance <- ggplot(sig_features, aes(x = OTU, y = EffectSize_r, fill = ColorGroup)) + 
	  geom_bar(stat = "identity") +
	  #geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
	  coord_flip() + 
	  theme_bw() + 
	  ylab("EffectSize r") +
	 ggtitle(paste(week, "Feature Importance")) +
	  scale_fill_manual(values = my_colors) #+ theme_classic()
	pImportance	
	
	saveRDS(file=paste0(week,".plot.rds"),pImportance)
	saveRDS(file=paste0(week,".sig_features.rds"),sig_features)
	write.table(file= paste0(week,".sig_features.tsv"),sig_features,row.names = F,sep="\t")
	return(pImportance)
}

###########################################################################################################################
##########################################             READ DATA                     ######################################
###########################################################################################################################

###############################
####  Read the infiles   ######

load("WT_w3.RData")
load("WT_w12.RData")


plot_w12 <- WTdiff_filter_plot( week = "w12", sig_features = WT_w12 )
plot_w3 <- WTdiff_filter_plot(  week = "w3", sig_features = WT_w3 )

ggsave( file = "wt_w3.pdf" , plot_w3, width = 8, height = 27 )
ggsave( file = "wt_w12.pdf" , plot_w12, width = 8, height = 27 )


