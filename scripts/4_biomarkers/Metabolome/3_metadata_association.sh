### sh 3_metadata_association.sh
metadata_association=/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/Metabolome/metadata_association_Metabolome.R

mkdir "Metadata_association"
cd Metadata_association
normfile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/Metabolomics/physeq.metabolites.norm.RData"
Variable="Disease"
SigF="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/Metabolome/DA/Disease/q_values_table.tsv"

Rscript --vanilla  $metadata_association $normfile $Variable $SigF
