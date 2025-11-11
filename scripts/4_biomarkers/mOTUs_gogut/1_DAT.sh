### sh 1_DAT.sh
DAT=./Differential_abundance_test_QMP.R

mkdir "DA"
cd DA

rarfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/Treatment_response/physeq.TreatRes.qmp.motus.RData"
Variable="Response"

Rscript --vanilla  $DAT $rarfile $Variable



