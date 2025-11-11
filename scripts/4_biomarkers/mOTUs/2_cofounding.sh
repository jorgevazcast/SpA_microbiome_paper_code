### sh 2_cofounding.sh
cofounding=./cofounding_QMP.R

mkdir "Cofounding"
cd Cofounding

rarfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/QMP_mOTUS/physeq.qmp.motus.RData"
Variable="Disease"
SigF="~/github_shared_code_and_publications/SpA_microbiome_paper_code/4_biomarkers/mOTUs/DA/Disease/q_values_table.tsv"

Rscript --vanilla  $cofounding $rarfile $Variable $SigF



