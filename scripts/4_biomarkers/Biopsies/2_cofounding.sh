### sh 2_cofounding.sh
cofounding=./cofounding.R

mkdir "Cofounding"
cd Cofounding

#Tax="ID"  ### Consensus_sp  silva_v138_2 ID
Tax=$1  ### Consensus_sp  silva_v138_2 ID
mkdir $Tax
cd $Tax
WD=$(pwd)

#### Ileum  ####
Biopsie="Ileum"
mkdir $Biopsie
cd $Biopsie

GMPRfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/"$Biopsie"_biopsies/GTDB_r220/physeq_sv.rar.gmpr.rds"
rarfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/"$Biopsie"_biopsies/GTDB_r220/physeq_sv.rar.rds"
above10000file="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/"$Biopsie"_biopsies/GTDB_r220/physeq_sv.rds"
clrfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/"$Biopsie"_biopsies/GTDB_r220/physeq_sv.rar.clr.rds"
Variable="Disease"
SigF="~/github_shared_code_and_publications/SpA_microbiome_paper_code/4_biomarkers/Biopsies/DA/$Tax/$Biopsie/$Variable/q_values_table.tsv"

Rscript --vanilla  $cofounding $GMPRfile $rarfile $above10000file $clrfile $Variable $Tax $SigF
cd $WD

#### Colon  ####
Biopsie="Colon"
mkdir $Biopsie
cd $Biopsie

GMPRfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/"$Biopsie"_biopsies/GTDB_r220/physeq_sv.rar.gmpr.rds"
rarfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/"$Biopsie"_biopsies/GTDB_r220/physeq_sv.rar.rds"
above10000file="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/"$Biopsie"_biopsies/GTDB_r220/physeq_sv.rds"
clrfile="~/github_shared_code_and_publications/SpA_microbiome_paper_code/1_infiles/"$Biopsie"_biopsies/GTDB_r220/physeq_sv.rar.clr.rds"
Variable="Disease"
SigF="~/github_shared_code_and_publications/SpA_microbiome_paper_code/4_biomarkers/Biopsies/DA/$Tax/$Biopsie/$Variable/q_values_table.tsv"

Rscript --vanilla  $cofounding $GMPRfile $rarfile $above10000file $clrfile $Variable $Tax $SigF
cd $WD


