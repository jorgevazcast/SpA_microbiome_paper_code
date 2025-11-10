### sh 1_DAT.sh
DAT=/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/mOTUs_gogut/Differential_abundance_test_QMP.R

mkdir "DA"
cd DA

rarfile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/Treatment_response/physeq.TreatRes.qmp.motus.RData"
Variable="Response"

Rscript --vanilla  $DAT $rarfile $Variable



