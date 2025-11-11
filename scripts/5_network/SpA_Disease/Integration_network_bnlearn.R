set.seed(12345)
library(phyloseq)
library("ggplot2")
library(ggpubr)
library(gridExtra)
library(microbiome)
library(igraph)
library(bnlearn)
library(glmmTMB)

path_func <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code/functions"

source(paste0( path_func,"/Network_supplementary_functions.R" ))
source(paste0( path_func,"/bnlearn_supplementary_functions.R" ))

#####################################################################################################################################################
##########################################################         SCRIPT             ###############################################################
#####################################################################################################################################################
path_wd <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles"

######################################
###### Read the phyloseq object ######

##  physeq.qmp.motus
load(paste0(path_wd,"/QMP_mOTUS/physeq.qmp.motus.RData"))
taxa_names(physeq.qmp.motus) <- gsub("[^[:alnum:]]", "_", taxa_names(physeq.qmp.motus))
taxa_names(physeq.qmp.motus) <- gsub(" ", "_", taxa_names(physeq.qmp.motus))	

##  physeq.qmp.kegg
#load("/home/luna.kuleuven.be/u0141268/Dropbox/1-GIANT_last/final_analysis/abundances/physeq.qmp.kegg.RData")

##  physeq.qmp.gmm
load(paste0(path_wd,"/QMP_GMM/physeq.qmp.gmm.RData"))
taxa_names(physeq.qmp.gmm) <- gsub("[^[:alnum:]]", "_", taxa_names(physeq.qmp.gmm))
taxa_names(physeq.qmp.gmm) <- gsub(" ", "_", taxa_names(physeq.qmp.gmm))	

##  physeq.metabolites.norm
load(paste0(path_wd,"/Metabolomics/physeq.metabolites.norm.RData"))
taxa_names(physeq.metabolites.norm) <- gsub("[^[:alnum:]]", "_", taxa_names(physeq.metabolites.norm))
taxa_names(physeq.metabolites.norm) <- gsub(" ", "_", taxa_names(physeq.metabolites.norm))	

#################################
######## Read biomarkers ########
biomarkers_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code/4_biomarkers"
infile <- paste0(biomarkers_dir,"/GMM/Metadata_association/Disease/rar_cofound_test_anova.tsv")
GMM.biomarkers <- read.table(infile,sep="\t",header = T)

infile <- paste0(biomarkers_dir,"/Metabolome/Metadata_association/Disease/log_cofound_test_anova.tsv")
Metabolome.biomarkers <- read.table(infile,sep="\t",header = T)

infile <- paste0(biomarkers_dir,"/mOTUs/Metadata_association/Disease/rar_cofound_test_anova.tsv")
mOTUs.biomarkers <- read.table(infile,sep="\t",header = T)

###################################
########  Filter samples   ########

Samples2Use <- table(c( sample_names(physeq.qmp.motus)  , sample_names(physeq.qmp.gmm), sample_names(physeq.metabolites.norm)  ))
Samples2Use <- names(Samples2Use[Samples2Use == 3])

physeq.qmp.motus <- prune_samples( Samples2Use , physeq.qmp.motus )
physeq.qmp.gmm <- prune_samples( Samples2Use , physeq.qmp.gmm )
physeq.metabolites.norm <- prune_samples( Samples2Use , physeq.metabolites.norm )

###################################
######## Filter biomarkers ########

GMM.biomarkers <- subset(GMM.biomarkers,q.value <= 0.1)
Metabolome.biomarkers <- subset(Metabolome.biomarkers,q.value <= 0.1)
mOTUs.biomarkers <- subset(mOTUs.biomarkers,q.value <= 0.1)

motus_df <- data.frame( otu_table(prune_taxa( unique(mOTUs.biomarkers$Feature) , physeq.qmp.motus )) )
gmm_df <- data.frame( otu_table(prune_taxa( unique(GMM.biomarkers$Feature) , physeq.qmp.gmm )) )
metabolites_df <- data.frame( otu_table(prune_taxa( unique(Metabolome.biomarkers$Feature) , physeq.metabolites.norm )) )

##################################################
######## Order the biomarkers data.frames ########
motus_df <- motus_df[,match( Samples2Use, colnames(motus_df))]
gmm_df <- gmm_df[,match( Samples2Use, colnames(gmm_df))]
metabolites_df <- metabolites_df[,match( Samples2Use, colnames(metabolites_df))]

##################################################
################     Metadata     ################
Sample_data <- data.frame(sample_data(physeq.qmp.motus))
Sample_data <- Sample_data[match( Samples2Use, rownames(Sample_data)),]

SigVar <- unique(c(GMM.biomarkers$Variable2Test,Metabolome.biomarkers$Variable2Test, mOTUs.biomarkers$Variable2Test))


Sample_data <- Sample_data[,match( SigVar, colnames(Sample_data) )]
Sample_data <- Sample_data[match( Samples2Use, rownames(Sample_data)),]


Disease <- binary_matrix(var=Sample_data$Disease,var.name="Disease",return_factor=T)
Disease.SpA <- Disease[,c("Disease.SpA")]
Disease_activity <- binary_matrix(var=Sample_data$Disease_activity,var.name="Disease_activity",return_factor=T)
Disease_activity <- Disease_activity[,c("Disease_activity.SpAHigh","Disease_activity.SpALow")]

NSAID_use <- binary_matrix(var=Sample_data$NSAID_use,var.name="NSAID_use",return_factor=T)
NSAID_use.NO <- NSAID_use[,c("NSAID_use.NO")]
Smoking <- binary_matrix(var=Sample_data$Smoking,var.name="Smoking",return_factor=T)
Smoking.history <- Smoking[,c("Smoking.history")]
inflam_back_pain <- binary_matrix(var=Sample_data$inflam_back_pain...47,var.name="inflam_back_pain...47",return_factor=T)
enthesitis_general <- binary_matrix(var=Sample_data$enthesitis_general,var.name="enthesitis_general",return_factor=T)
uveitis <- binary_matrix(var=Sample_data$uveitis,var.name="uveitis",return_factor=T)

Fecal_calpro_values <- Sample_data$Fecal_calpro_values
Age_at_visit <- Sample_data$Age_at_visit

##########################################################
################     Join the tables      ################


QMP <- rbind(motus_df,gmm_df)
QMP <- t(QMP)
colnames(QMP)  <- paste0("QMP.",colnames(QMP))


Linear  <- rbind(Age_at_visit=Age_at_visit,Fecal_calpro_values=Fecal_calpro_values,metabolites_df)
Linear <- t(Linear)
colnames(Linear)  <- paste0("Linear.",colnames(Linear))


Boolean  <- data.frame(Disease.SpA=Disease.SpA,
		Disease_activity,NSAID_use.NO=NSAID_use.NO,
		Smoking.history=Smoking.history, enthesitis_general,uveitis,inflam_back_pain)
colnames(Boolean)  <- paste0("Boolean.",colnames(Boolean))
rownames(Boolean) <- rownames(Sample_data)


df2bnlearn <-cbind(Boolean,Linear,QMP)
df2bnlearn<-df2bnlearn[complete.cases(df2bnlearn),]

retHC <- hc( df2bnlearn, score = "custom", fun = my.bic.GLM, args = list()   )
retBNlearns <- retHC



##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################

