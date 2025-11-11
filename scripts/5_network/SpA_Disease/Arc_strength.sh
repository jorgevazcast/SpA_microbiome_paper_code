#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pe serial 111  
#$ -l mem_limit=200G

source /etc/profile.d/modules.sh
module load R/4.0.3

Rscript --vanilla Arc_strength.R

# qsub Arc_strength.R

