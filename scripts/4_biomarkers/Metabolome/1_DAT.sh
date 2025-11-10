### sh 1_DAT.sh
DAT=/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/Metabolome/Differential_abundance_test_Metabolome.R

mkdir "DA"
cd DA
#Tax="ID"  ### Consensus_sp  silva_v138_2 ID
normfile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/Metabolomics/physeq.metabolites.norm.RData"
rawfile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/Metabolomics/physeq.metabolites.RData"
Variable="Disease"

Rscript --vanilla  $DAT $normfile $rawfile $Variable



