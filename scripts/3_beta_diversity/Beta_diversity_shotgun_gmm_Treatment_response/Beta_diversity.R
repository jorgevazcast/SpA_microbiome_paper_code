# Rscript --vanilla Beta_diversity.R
set.seed(12345)
library(phyloseq)
library("ggplot2")
library(ggpubr)
library(gridExtra)

working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 
source(paste0(path_dir, "/functions/beta_diver_functions_V9.R"))
source(paste0(path_dir, "/functions/supplementary_functions_beta_diver.R"))

#####################################################################################################################################################
##########################################################       SCRIPT         #####################################################################
#####################################################################################################################################################

FDR_pval <- 0.1

######################################
###### Read the phyloseq object ######
load(paste0(path_dir, "/1_infiles/Treatment_response/physeq.TreatRes.qmp.gmm.RData"))

physeq_QMP <- physeq.TreatRes.qmp.gmm
otable <- as.matrix( otu_table(physeq_QMP))

#####################################
###### Filter the taxa (motus) ######
taxa_QMP<-filter_low_prevalence(in.table = otable, max_percentage_0 = 80 ) 

###############################
###### Variables to use  ######
var2se <- c("ASDAS_W0","BASDAI_tot_W0","Gut_inflammation","Fecal_calpro_values","Treatment_response","Colon_inflammation","Ileum_inflammation")

####### Subset the variables #######
Metadata  <- data.frame( sample_data(physeq_QMP) , stringsAsFactors = F )
Metadata <- Metadata[,match(var2se,colnames(Metadata))]

PercentNA <- sapply(colnames(Metadata), function(x,Met2=Metadata){
		NA_N <- table(is.na(Met2[,x]))
		if( any(names(NA_N) == "TRUE") ){
			print(x)
			ret <- NA_N["TRUE"]/nrow(Met2) * 100
		} else{ret <- 0}
		return(ret)	
	}
)

####### Subset the variables by the percentange of NA #######
names(PercentNA) <- gsub(".TRUE","",names(PercentNA))
cbind( colnames(Metadata)  , names(PercentNA)    )
Varaibles2use <- names(PercentNA[PercentNA <= 10])

####### Check the varaible class #######
for(i in Varaibles2use){
	print(i)
	print(class(Metadata[,i]))
	print(length(unique(Metadata[,i])))
	
	if(length(unique(Metadata[,i])) > 5 &  class(Metadata[,i]) == "character"){
		print(i)
		Metadata[,i] <- as.numeric(as.character(Metadata[,i]))
	}
}

#################################################################################
########################	ADONIS		###############################
#################################################################################
Metadata2use <- Metadata[,match(Varaibles2use,colnames(Metadata))]
Metadata2use$Response <- ifelse(Metadata2use$Treatment_response == "Non_responder" , "Non_responder","Responder")
Metadata2use$Treatment_response <- NULL
taxa_QMP <- taxa_QMP[,match(rownames(Metadata2use), colnames(taxa_QMP))]
taxa_QMP <- taxa_QMP[rowSums(taxa_QMP) != 0,]
taxa_matrix <- t(taxa_QMP)

out_dir  = "./"
bray.table.ADONIS <- ADONIS_func( in.matrix =  taxa_matrix  , Distance = "bray", in.Metadata = data.frame(Metadata2use) , prefix = paste0(out_dir,"/bray") )
bray.table.ADONIS <- bray.table.ADONIS[order(bray.table.ADONIS$p.value),]


#################################################################################
########################	Plot all the PCoA	###############################
#################################################################################
PCoA_plot_dir <- "./PCoA_plot_dir"
dir.create(PCoA_plot_dir)
list_PCoA <- list()
sum <- 1
for(i in rownames(bray.table.ADONIS)){
	PCoA <- PCoA_grapper(Var = i, tempMetadata = Metadata2use, bray.table.ADONIS=bray.table.ADONIS,taxa_mat=taxa_matrix)
	list_PCoA[[sum]] <-PCoA
	ggsave(paste0(PCoA_plot_dir,"/",i,".pdf"), PCoA, width=7, height=6)
	sum <- sum +1
	rm(PCoA	)


}

AllPCoAs  <-   ggarrange( plotlist = list_PCoA )
ggsave("AllPCoAs.pdf", AllPCoAs, width=40, height=20)
AllPCoAs

#################################################################################
########################	Plot all the NMDS	###############################
#################################################################################
NMDS_plot_dir <- "./NMDS_plot_dir"
dir.create(NMDS_plot_dir)
list_NMDS <- list()
sum <- 1
for(i in rownames(bray.table.ADONIS)){
	NMDS <- nmds_grapper(Var = i, tempMetadata = Metadata2use, bray.table.ADONIS=bray.table.ADONIS,taxa_mat=taxa_matrix)
	list_NMDS[[sum]] <-NMDS
	ggsave(paste0(NMDS_plot_dir,"/",i,".pdf"), NMDS, width=7, height=6)
	sum <- sum +1
	rm(NMDS	)


}

AllNMDS  <-   ggarrange( plotlist = list_NMDS )
ggsave("AllNMDS.pdf", AllNMDS, width=40, height=20)
AllNMDS


#################################################################################
########################	ordiR2step		#####################
#################################################################################
bray.table.ADONIS <- data.frame(bray.table.ADONIS)
bray.table.ADONIS$BH.adj.p.value <- as.numeric(as.character(bray.table.ADONIS$BH.adj.p.value))
bray.table.ADONIS <- subset(bray.table.ADONIS,BH.adj.p.value < FDR_pval)

Metadata2use <- Metadata2use[ , match( as.character(bray.table.ADONIS$Variable) , colnames(Metadata2use)) ]
Metadata2use <- Metadata2use[complete.cases(Metadata2use),]
N<-dim(Metadata2use)[1] # N samples

#####  bray
taxa_matrix <- taxa_matrix[match( rownames(Metadata2use) ,  rownames(taxa_matrix)),]

########################	ordiR2step		#########################
capscale_bray<-capscale_cum_variance( in.matrix =  taxa_matrix , Distance = "bray", in.Metadata = data.frame(Metadata2use),adj.pval.cutof =FDR_pval, prefix = "bray")
capscale_bray
#capscale_bray[["non_redundant"]]

#### write tables 
Table_capscale_all <- data.frame(capscale_bray$all)
Table_capscale_all[order(Table_capscale_all$p.value),]
Table_capscale_all$N <- N
write.table(Table_capscale_all,"Table_capscale_all.tsv",col.names=T,row.names = T,quote=FALSE,sep = "\t")

Table_capscale_non_redundant <- data.frame(capscale_bray$non_redundant)
Table_capscale_non_redundant$N <- N
write.table(Table_capscale_non_redundant,"Table_capscale_non_redundant.tsv",col.names=T,row.names = T,quote=FALSE,sep = "\t")


#####  Plot  Bray 
capscale_bray2plot<-nonredundant2plot(non.redundant = capscale_bray[["non_redundant"]], all.var = capscale_bray[["all"]] )

P<- ggplot(data=capscale_bray2plot[["nr"]], aes(x=Variable, y=R2, fill=Variance)) +
	geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	 ggtitle("BC distance all variables")
ggsave("Bray_ordiR2step.pdf",plot=P)	


if(all(capscale_bray != "Non significant variables")){
	NR_capscale <- capscale_bray[["non_redundant"]]
	capscale_bray2plot<-nonredundant2plot(non.redundant = NR_capscale, all.var = capscale_bray[["all"]] )
	capscale_bray2plot_sig <- capscale_bray2plot[["nr_sig"]]

	Pclr<- ggplot(data=capscale_bray2plot_sig, aes(x=Variable, y=R2, fill=Variance)) + geom_bar(stat="identity", position=position_dodge())  +
		theme_minimal() + ggtitle("Non-redundant variation") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))    #+  coord_flip() 
	ggsave("Bray_ordiR2step_sig.pdf",plot=Pclr,width = 12)	

	capscale_bray2plot_sig <- capscale_bray2plot_sig[capscale_bray2plot_sig$Variance != "Unconstrain",]
	capscale_bray2plot_sig$Variance <- factor(as.character(capscale_bray2plot_sig$Variance), levels=c("R2","Cum_R2"))
	PvarNR <- ggplot(data=capscale_bray2plot_sig, aes(x=Variable, y=R2, fill=Variance)) + geom_bar(stat="identity", position=position_dodge())  +
		theme_minimal() + ggtitle("Non-redundant variation") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) + #+  coord_flip() 
		scale_fill_manual(values=c("#440154","#21908c"))
	ggsave("Bray_ordiR2step_sig_brief.pdf",plot=PvarNR,width = 7)
		
	
}


