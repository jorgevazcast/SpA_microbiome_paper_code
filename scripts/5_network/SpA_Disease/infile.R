set.seed(12345)
library(phyloseq)
library(microbiome)
library(igraph)
library(bnlearn)

path_func <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code/functions"
source(paste0( path_func,"/Functions/Extended_Statistical_Toolkit_functions.R" ))
source(paste0( path_func,"/Network_supplementary_functions.R" ))
source(paste0( path_func,"/bnlearn_supplementary_functions.R" ))

#####################################################################################################################################################
##########################################################         SCRIPT             ###############################################################
#####################################################################################################################################################
path_wd <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles"
Nboot <- 1000
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
biomarkers_dir <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/Cofound_glmnet_biomarkers"
infile <- paste0(biomarkers_dir,"/Log/MeanImp.RData")
load(infile)

MeanImp <- MeanImp[MeanImp$Freq > round(Nboot * 0.7),]


###################################
########  Filter samples   ########

Samples2Use <- table(c( sample_names(physeq.qmp.motus)  , sample_names(physeq.qmp.gmm), sample_names(physeq.metabolites.norm)  ))
Samples2Use <- names(Samples2Use[Samples2Use == 3])

physeq.qmp.motus <- prune_samples( Samples2Use , physeq.qmp.motus )
physeq.qmp.gmm <- prune_samples( Samples2Use , physeq.qmp.gmm )
physeq.metabolites.norm <- prune_samples( Samples2Use , physeq.metabolites.norm )

###################################
######## Filter biomarkers ########

mOTUsbio <- taxa_names(physeq.qmp.motus)
mOTUsbio <- mOTUsbio[mOTUsbio %in% as.character(MeanImp$Feature)]

gmmbio <- taxa_names(physeq.qmp.gmm)
gmmbio <- gmmbio[gmmbio %in% as.character(MeanImp$Feature)]

metbio <- taxa_names(physeq.metabolites.norm)
metbio <- metbio[metbio %in% as.character(MeanImp$Feature)]


motus_df <- data.frame( otu_table(prune_taxa( mOTUsbio , physeq.qmp.motus )) )
gmm_df <- data.frame( otu_table(prune_taxa( gmmbio , physeq.qmp.gmm )) )
metabolites_df <- data.frame( otu_table(prune_taxa( metbio , physeq.metabolites.norm )) )

##################################################
######## Order the biomarkers data.frames ########
motus_df <- motus_df[,match( Samples2Use, colnames(motus_df))]
gmm_df <- gmm_df[,match( Samples2Use, colnames(gmm_df))]
metabolites_df <- metabolites_df[,match( Samples2Use, colnames(metabolites_df))]

##################################################
################     Metadata     ################
Sample_data <- data.frame(sample_data(physeq.qmp.motus))
Sample_data <- Sample_data[match( Samples2Use, rownames(Sample_data)),]

SigVar <- colnames(Sample_data)
SigVar <- SigVar[SigVar %in% as.character(MeanImp$Feature)]
if( any(  as.character(MeanImp$Feature) %in% "Bacteroides_2" )){
	SigVar <- unique(c(SigVar,"enterotype"))
}
SigVar <- unique(c("Disease",SigVar))

Sample_data <- Sample_data[,match( SigVar, colnames(Sample_data) )]
Sample_data <- Sample_data[match( Samples2Use, rownames(Sample_data)),]

### Boolean
Disease <- binary_matrix(var=Sample_data$Disease,var.name="Disease",return_factor=T)
Disease.SpA <- Disease[,c("Disease.SpA")]

Enterotype <- binary_matrix(var=Sample_data$enterotype,var.name="enterotype",return_factor=T)
##########################################################
################     Join the tables      ################


QMP <- rbind(motus_df,gmm_df)
QMP <- t(QMP)
colnames(QMP)  <- paste0("QMP.",colnames(QMP))


Linear  <- rbind(
#		Age_at_visit=Age_at_visit,
#		Fecal_calpro_values=Fecal_calpro_values,
#		Water=Water,
#		BMI=BMI,
		metabolites_df)
Linear <- t(Linear)
colnames(Linear)  <- paste0("Linear.",colnames(Linear))


Boolean  <- data.frame(
		Disease.SpA=Disease.SpA,
#		Disease_activity,
#		NSAID_use.NO=NSAID_use.NO,
#		Smoking.history=Smoking.history, 
#		enthesitis_general,
#		uveitis,inflam_back_pain,
#		Enterotype)
		Bacteroides_2 = Enterotype$enterotype.Bacteroides_2 
	)
		
colnames(Boolean)  <- paste0("Boolean.",colnames(Boolean))
rownames(Boolean) <- rownames(Sample_data)


df2bnlearn <-cbind(Boolean,Linear,QMP)


#######################################################################
######## Summarizes the quality and structure of each variable ########
#######################################################################

#### Remove the varaibles that only have a single category or with a high number of NAs ####
N_values_NA_cat <- N_categories_NA_values_function( in.meta = df2bnlearn )

### High NA data
High_NA_var <- N_values_NA_cat[N_values_NA_cat$NA_percentage >= 20,]$Variable

### Single categorie variable
Single_categorie <- subset(N_values_NA_cat,  N_values == 1 )$Variable

### Remove binary variables with less than 1% of imbalance data
min_num_samp_per_cat <- round(nrow(df2bnlearn) * 0.1)
Non_Info_var <- subset(N_values_NA_cat, Var_type == "factor" & N_min_categories < min_num_samp_per_cat & N_values == 2 )$Variable

#### Exclude the varaibles ####
Exclude_variables <- unique( c(    High_NA_var, Single_categorie, Non_Info_var ) )
df2bnlearn <- df2bnlearn[,!colnames(df2bnlearn) %in% Exclude_variables]
dim(df2bnlearn)

# df2bnlearn <-df2bnlearn[,c("Boolean.Disease.SpA",subVars)]
df2bnlearn<-df2bnlearn[complete.cases(df2bnlearn),]


#######################################################################
######## Summarizes the quality and structure of each variable ########
#######################################################################

dir.create("infile")
saveRDS(file = "./infile/df2bnlearn.rds" , df2bnlearn)

