#!/bin/bash
#$ -t 1-500
#$ -cwd
#$ -S /bin/bash
##$ -pe serial 1  
#$ -l mem_limit=3G

source /etc/profile.d/modules.sh
module load R/4.0.3

Rscript --vanilla Boot_bnlearn.R $SGE_TASK_ID

# qsub Boot_bnlearn.sh

