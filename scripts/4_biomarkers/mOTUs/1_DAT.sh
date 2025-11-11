### sh 1_DAT.sh
DAT=./Differential_abundance_test_QMP.R

mkdir "DA"
cd DA
#Tax="ID"  ### Consensus_sp  silva_v138_2 ID
rarfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/QMP_mOTUS/physeq.qmp.motus.RData"
Variable="Disease"

Rscript --vanilla  $DAT $rarfile $Variable



