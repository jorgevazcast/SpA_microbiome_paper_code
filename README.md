# SpA_microbiome_paper_code
This repository contains a comprehensive pipeline for analyzing gut microbiome data from the Spondyloarthritis (SpA) cohort. The workflow includes all major steps of the analysis, including data preprocessing, quality control, taxonomic and functional profiling, statistical modeling, and visualization. The documentation below provides a detailed guide on how to reproduce each step of the analysis.

## Software requirements
###### R version 4.1.2 (2021-11-01) -- "Bird Hippie"
R packages:
* caret_6.0-86
* DMwR_0.4.1
* ROSE_0.0-4
* pROC_1.17.0.1
* mltools_0.3.5

Important! Change the working directory before loading the different functions
```R
working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code"
path_func <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code/functions" 
```

## Scripts
All scripts used in the analysis are available in the `script` directory
>
**2_alpha_diversity**
* `/scripts/2_alpha_diversity/Colon/Diversity_test.R`
* `/scripts/2_alpha_diversity/GMM/Diversity_test.R`
* `/scripts/2_alpha_diversity/Ileum/Diversity_test.R`
* `/scripts/2_alpha_diversity/mOTUs/Diversity_test.R`

**3_beta_diversity**
* `/scripts/3_beta_diversity/Beta_diversity_colon_biopsies_genus/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_colon_biopsies_sv/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_ileum_biopsies_sv/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_metabolomics/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_shotgun_gmm_Treatment_response/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_shotgun_motus/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_shotgun_motus_Treatment_response/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_shotgun_gmm/Bray/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_shotgun_gmm/Canberra/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_shotgun_kegg/Bray/Beta_diversity.R`
* `/scripts/3_beta_diversity/Beta_diversity_shotgun_kegg/Canberra/Beta_diversity.R`

**4_biomarkers**
>
Biomarkers
* `/scripts/4_biomarkers/Biomarkers_response_gogut_prediction/boot632_biomarkers.R`
* `/scripts/4_biomarkers/Biomarkers_response_gogut_stability/glmnet.R`
* `/scripts/4_biomarkers/Cofound_glmnet_biomarkers/glmnet.R`
>
Differential abundance
###### Run the Differential abundance in the BEGiant dataset
* `/scripts/4_biomarkers/Biopsies/Biopsies_DA_pipeline.sh`
* `/scripts/4_biomarkers/GMM/Biopsies_DA_pipeline.sh`
* `/scripts/4_biomarkers/mOTUs/Biopsies_DA_pipeline.sh`
>
###### Differential abundance in the gogut dataset
* `/scripts/4_biomarkers/Metabolome/Biopsies_DA_pipeline.sh`
* `/scripts/4_biomarkers/GMM_gogut/1_DAT.sh`
* `/scripts/4_biomarkers/mOTUs_gogut/1_DAT.sh`

**5_Bayesian_network**
>
Bayesian network inference
* `/scripts/5_network/SpA_Disease/Bayesian_network_pipeline.sh`
>
Bootstrapping and arc strength
###### The arc strength estimation and Bayesian network bootstrapping were performed using a Sun Grid Engine (SGE) queuing cluster architecture (via a qsub submission script)
* `/scripts/5_network/SpA_Disease/Cluster_scripts/run_bn_bootstrap_learning.sh`

**6_mice_experiments**
>
/home/luna.kuleuven.be/u0141268/github_shared_code_and_publications/SpA_microbiome_paper_code/scripts/6_mice_experiments/Diff_abundance_Qfemto/run_Diff_abundance_Qfemto.sh


## Raw data
FASTQ files for the project can be found under the following EGA accession number



