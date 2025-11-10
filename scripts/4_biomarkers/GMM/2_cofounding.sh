### sh 2_cofounding.sh
cofounding=/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/GMM/cofounding_QMP.R

mkdir "Cofounding"
cd Cofounding

rarfile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/QMP_GMM/physeq.qmp.gmm.RData"
Variable="Disease"
SigF="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/GMM/DA/Disease/q_values_table.tsv"

Rscript --vanilla  $cofounding $rarfile $Variable $SigF



