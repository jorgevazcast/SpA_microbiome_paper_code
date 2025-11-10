
### DA ###
sh 1_DAT.sh Consensus_sp
sh 1_DAT.sh silva_v138_2
sh 1_DAT.sh ID

### cofounding ####
sh 2_cofounding.sh Consensus_sp
sh 2_cofounding.sh silva_v138_2
sh 2_cofounding.sh ID

#### Metadata association ####
sh 3_metadata_association.sh Consensus_sp
sh 3_metadata_association.sh silva_v138_2
sh 3_metadata_association.sh ID

