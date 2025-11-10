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
```

## Scripts
All scripts used in the analysis are available in the `script` directory
>
**2_alpha_diversity**
###### /scripts/2_alpha_diversity/Colon
###### `/scripts/2_alpha_diversity/GMM`
* `/scripts/2_alpha_diversity/Ileum`
* `/scripts/2_alpha_diversity/mOTUs`

## Raw data
FASTQ files for the project can be found under the following EGA accession number



