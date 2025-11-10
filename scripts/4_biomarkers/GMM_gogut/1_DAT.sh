### sh 1_DAT.sh
DAT=/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/GMM_gogut/Differential_abundance_test_QMP.R

mkdir "DA"
cd DA
#Tax="ID"  ### Consensus_sp  silva_v138_2 ID
rarfile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/Treatment_response/physeq.TreatRes.qmp.gmm.RData"
Variable="Response"

Rscript --vanilla  $DAT $rarfile $Variable



