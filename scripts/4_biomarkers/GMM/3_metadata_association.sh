### sh 3_metadata_association.sh
metadata_association=/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/GMM/metadata_association_QMP.R

mkdir "Metadata_association"
cd Metadata_association
rarfile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/QMP_GMM/physeq.qmp.gmm.RData"
Variable="Disease"
SigF="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/GMM/DA/Disease/q_values_table.tsv"

Rscript --vanilla  $metadata_association $rarfile $Variable $SigF
