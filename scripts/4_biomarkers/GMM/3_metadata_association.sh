### sh 3_metadata_association.sh
metadata_association=./metadata_association_QMP.R

mkdir "Metadata_association"
cd Metadata_association
rarfile="~/github_shared_code_and_publications/1_infiles/QMP_GMM/physeq.qmp.gmm.RData"
Variable="Disease"
SigF="~/github_shared_code_and_publications/4_biomarkers/GMM/DA/Disease/q_values_table.tsv"

Rscript --vanilla  $metadata_association $rarfile $Variable $SigF
