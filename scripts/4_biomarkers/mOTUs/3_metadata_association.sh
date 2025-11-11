### sh 3_metadata_association.sh
metadata_association=./metadata_association_QMP.R

mkdir "Metadata_association"
cd Metadata_association
rarfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/QMP_mOTUS/physeq.qmp.motus.RData"
Variable="Disease"
SigF="~/github_shared_code_and_publications/SpA_microbiome_paper_code/4_biomarkers/mOTUs/DA/Disease/q_values_table.tsv"

Rscript --vanilla  $metadata_association $rarfile $Variable $SigF
