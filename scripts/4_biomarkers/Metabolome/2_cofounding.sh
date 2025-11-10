### sh 2_cofounding.sh
cofounding=/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/Metabolome/cofounding_Metabolome.R

mkdir "Cofounding"
cd Cofounding

normfile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/Metabolomics/physeq.metabolites.norm.RData"
Variable="Disease"
SigF="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/Metabolome/DA/Disease/q_values_table.tsv"

Rscript --vanilla  $cofounding $normfile $Variable $SigF



