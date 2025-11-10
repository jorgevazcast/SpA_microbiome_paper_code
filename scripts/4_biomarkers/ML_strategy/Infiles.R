set.seed(12345)
library(phyloseq)
library(microbiome)

path_func <- "/home/luna.kuleuven.be/u0141268/github_projects/phyloseq_format_converter"
source(paste0(path_func,"/Functions/Utils.R"))


#####################################################################################################################################################
##########################################################         SCRIPT             ###############################################################
#####################################################################################################################################################
path_wd <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles"

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
biomarkers_dir <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers"
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


##################################################
################     Outfiles     ################
in_data <- rbind( motus_df  , gmm_df  ,metabolites_df)
log_data <- rbind(log(motus_df + 1),log(gmm_df + 1),metabolites_df)
scaled_data <- rbind( t(scale(t(motus_df)))  ,  t(scale(t(gmm_df)))  ,  metabolites_df)

FakeTax <- as.matrix(data.frame( Dummy = "N", Taxa = rownames(scaled_data)    ))
rownames(FakeTax)  <- rownames(scaled_data)

phylo_biomarkers_log <- create_phyloseq_object_fun( data = log_data, tax = FakeTax, sample_dat = Sample_data  )
phylo_biomarkers_scaled <- create_phyloseq_object_fun( data = scaled_data, tax = FakeTax, sample_dat = Sample_data  )
phylo_biomarkers <- create_phyloseq_object_fun( data = in_data, tax = FakeTax, sample_dat = Sample_data  )

dir.create("infiles")
saveRDS(file="./infiles/phylo_biomarkers_log.rds",phylo_biomarkers_log)
saveRDS(file="./infiles/phylo_biomarkers.rds",phylo_biomarkers)
saveRDS(file="./infiles/phylo_biomarkers_scaled.rds",phylo_biomarkers_scaled)

#sex 
#age
#AST
#ALT
#category
