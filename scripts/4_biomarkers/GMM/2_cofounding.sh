### sh 2_cofounding.sh
cofounding=./cofounding_QMP.R

mkdir "Cofounding"
cd Cofounding

rarfile="~/github_shared_code_and_publications/1_infiles/QMP_GMM/physeq.qmp.gmm.RData"
Variable="Disease"
SigF="~/github_shared_code_and_publications/4_biomarkers/GMM/DA/Disease/q_values_table.tsv"

Rscript --vanilla  $cofounding $rarfile $Variable $SigF



