# Rscript --vanilla Beta_diversity.R
set.seed(12345)
library(phyloseq)
library("ggplot2")
library(ggpubr)
library(gridExtra)

working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

source(paste0(working_dir, "/functions/beta_diver_functions_V9.R"))
source(paste0(working_dir, "/functions/supplementary_functions_beta_diver.R"))

#####################################################################################################################################################
##########################################################       SCRIPT         #####################################################################
#####################################################################################################################################################

FDR_pval <- 0.1

######################################
###### Read the phyloseq object ######
load(paste0(working_dir, "/1_infiles/Metabolomics/physeq.metabolites.norm.RData"))


physeq_data <- physeq.metabolites.norm
otable <- as.matrix( otu_table(physeq_data))

#####################################
###### Filter the taxa (motus) ######
taxa_data<-filter_low_prevalence(in.table = otable, max_percentage_0 = 80 ) 

###############################
###### Variables to use  ######
var2se <- as.character(Variables2Use$Var)

load(paste0(working_dir, "/Metadata/metagenomic_varaibles_associations/Shotgun_varaibles/Var2eliminate.RData" ))

Variables2eliminate <- Var2eliminate$Var
Variables2eliminate <- c(Variables2eliminate,c("uveitis","enthesitis_general"))
var2se <- var2se[!var2se %in% Variables2eliminate]

####### Subset the variables #######
Metadata  <- data.frame( sample_data(physeq_data) , stringsAsFactors = F )
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
}

#################################################################################
########################	ADONIS		###############################
#################################################################################
Metadata2use <- Metadata[,match(Varaibles2use,colnames(Metadata))]
taxa_data <- taxa_data[,match(rownames(Metadata2use), colnames(taxa_data))]
taxa_data <- taxa_data[rowSums(taxa_data) != 0,]
taxa_matrix <- t(taxa_data)

out_dir  = "./"
euclidean.table.ADONIS <- ADONIS_func( in.matrix =  taxa_matrix  , Distance = "euclidean", in.Metadata = data.frame(Metadata2use) , prefix = paste0(out_dir,"/euclidean") )
euclidean.table.ADONIS <- euclidean.table.ADONIS[order(euclidean.table.ADONIS$p.value),]


#################################################################################
########################	Plot all the PCoA	###############################
#################################################################################
PCoA_plot_dir <- "./PCoA_plot_dir"
dir.create(PCoA_plot_dir)
list_PCoA <- list()
sum <- 1
for(i in rownames(euclidean.table.ADONIS)){
	PCoA <- PCoA_grapper(Var = i, tempMetadata = Metadata2use, bray.table.ADONIS=euclidean.table.ADONIS,taxa_mat=taxa_matrix,Distance="euclidean")
	list_PCoA[[sum]] <-PCoA
	ggsave(paste0(PCoA_plot_dir,"/",i,".pdf"), PCoA, width=7, height=6)
	sum <- sum +1
	rm(PCoA	)


}

AllPCoAs  <-   ggarrange( plotlist = list_PCoA )
ggsave("AllPCoAs.pdf", AllPCoAs, width=40, height=20)
AllPCoAs

#################################################################################
########################	ordiR2step		#####################
#################################################################################
euclidean.table.ADONIS <- data.frame(euclidean.table.ADONIS)
euclidean.table.ADONIS$BH.adj.p.value <- as.numeric(as.character(euclidean.table.ADONIS$BH.adj.p.value))
euclidean.table.ADONIS <- subset(euclidean.table.ADONIS,BH.adj.p.value < FDR_pval)

Metadata2use <- Metadata2use[ , match( as.character(euclidean.table.ADONIS$Variable) , colnames(Metadata2use)) ]
Metadata2use <- Metadata2use[complete.cases(Metadata2use),]
N<-dim(Metadata2use)[1] # N samples

#####  euclidean
taxa_matrix <- taxa_matrix[match( rownames(Metadata2use) ,  rownames(taxa_matrix)),]

########################	ordiR2step		#########################
capscale_euclidean<-capscale_cum_variance( in.matrix =  taxa_matrix , Distance = "euclidean", in.Metadata = data.frame(Metadata2use),adj.pval.cutof =FDR_pval, prefix = "euclidean")
capscale_euclidean
#capscale_euclidean[["non_redundant"]]

#### write tables 
Table_capscale_all <- data.frame(capscale_euclidean$all)
Table_capscale_all[order(Table_capscale_all$p.value),]
Table_capscale_all$N <- N
write.table(Table_capscale_all,"Table_capscale_all.tsv",col.names=T,row.names = T,quote=FALSE,sep = "\t")

Table_capscale_non_redundant <- data.frame(capscale_euclidean$non_redundant)
Table_capscale_non_redundant$N <- N
write.table(Table_capscale_non_redundant,"Table_capscale_non_redundant.tsv",col.names=T,row.names = T,quote=FALSE,sep = "\t")


#####  Plot  euclidean 
capscale_euclidean2plot<-nonredundant2plot(non.redundant = capscale_euclidean[["non_redundant"]], all.var = capscale_euclidean[["all"]] )

P<- ggplot(data=capscale_euclidean2plot[["nr"]], aes(x=Variable, y=R2, fill=Variance)) +
	geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	 ggtitle("BC distance all variables")
ggsave("euclidean_ordiR2step.pdf",plot=P)	


if(all(capscale_euclidean != "Non significant variables")){
	NR_capscale <- capscale_euclidean[["non_redundant"]]
	capscale_euclidean2plot<-nonredundant2plot(non.redundant = NR_capscale, all.var = capscale_euclidean[["all"]] )
	capscale_euclidean2plot_sig <- capscale_euclidean2plot[["nr_sig"]]

	Pclr<- ggplot(data=capscale_euclidean2plot_sig, aes(x=Variable, y=R2, fill=Variance)) + geom_bar(stat="identity", position=position_dodge())  +
		theme_minimal() + ggtitle("Non-redundant variation") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))    #+  coord_flip() 
	ggsave("euclidean_ordiR2step_sig.pdf",plot=Pclr,width = 12)	

	capscale_euclidean2plot_sig <- capscale_euclidean2plot_sig[capscale_euclidean2plot_sig$Variance != "Unconstrain",]
	capscale_euclidean2plot_sig$Variance <- factor(as.character(capscale_euclidean2plot_sig$Variance), levels=c("R2","Cum_R2"))
	PvarNR <- ggplot(data=capscale_euclidean2plot_sig, aes(x=Variable, y=R2, fill=Variance)) + geom_bar(stat="identity", position=position_dodge())  +
		theme_minimal() + ggtitle("Non-redundant variation") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) + #+  coord_flip() 
		scale_fill_manual(values=c("#440154","#21908c"))
	ggsave("euclidean_ordiR2step_sig_brief.pdf",plot=PvarNR,width = 7)
		
	
}


