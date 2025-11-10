#!/bin/bash
##$ -t 1-10
##$ -t 1
#$ -cwd
#$ -S /bin/bash
#$ -pe serial 110  
#$ -l mem_limit=70G
##$ -o ./log/log_$JOB_ID_$TASK_ID.out
##$ -e ./log/err_$TASK_ID.out

## qsub Classification_methods_cluster_Enterotype.sh
## sh Classification_methods.sh
#module load R/4.4.1
#source /etc/profile.d/modules.sh
#module load R/4.0.3

####################################################################################
#########################         Infile                   #########################
infile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/QMP_mOTUS/physeq.qmp.motus.RData"
infile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/QMP_GMM/physeq.qmp.gmm.RData"
infile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles/Metabolomics/physeq.metabolites.norm.RData"

infile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/ML_strategy/infiles/phylo_biomarkers.rds"
infile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/ML_strategy/infiles/phylo_biomarkers_log.rds"
infile="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/4_biomarkers/ML_strategy/infiles/phylo_biomarkers_scaled.rds"

######################################################################################################
#########################         Varaible, cohort, and working dir         #########################
COHORT="biomarkers" # mOTUs  GMM  Metabolome biomarkers
Variable="Disease"
Normalization="scaled" ### Rarefaction  CLR GMPR QMP Log

echo $Variable

### Set the directory
mkdir PREDICTIONS
cd PREDICTIONS
mkdir $COHORT
cd $COHORT
mkdir $Variable
cd $Variable
mkdir $Normalization
cd $Normalization

####################################################################################
#########################         Global variables         #########################
#HOMEC=/raeslab/scratch/jorvaz
HOMEC=$HOME

DIR_FUNCTIONS=$HOMEC"/github_projects/mlpredictr/Functions" ## MODIFY THIS VARIABLE WHERE YOU HAVE THE REPO

#### Train the models and choose the ones with better hyperparameters
ML_method_rf=$HOMEC"/github_projects/mlpredictr/Scripts/RandomForest.R"
ML_method_glmnet=$HOMEC"/github_projects/mlpredictr/Scripts/glmnet.R"
ML_method_XGboost=$HOMEC"/github_projects/mlpredictr/Scripts/XGboost.R"

### Script for obtaining the stats for each model
Plot_AUC_importance=$HOMEC"/github_projects/mlpredictr/Scripts/Plot_AUC_importance.R"

### Train the complete dataset using the hyperparametesr of the best model
TrainModel_Complete_dataset=$HOMEC"/github_projects/mlpredictr/Scripts/TrainModel_Complete_dataset.R"

### Prevalence
PrevCutoff=0.1

### Ncores
ncores=10



################################################################################
#########################            GLMNET            #########################
Feature_selection=F

### Train the hyperparameters ###
Rscript --vanilla $ML_method_glmnet $infile $Variable $PrevCutoff $Feature_selection $ncores $DIR_FUNCTIONS

### Stats best model ###
Dir_Results="./glmnet_FS_FALSE/"
Rscript --vanilla $Plot_AUC_importance $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS

### Train the best model using the best hyperparameters ###
Rscript --vanilla $TrainModel_Complete_dataset $infile $Variable $Dir_Results $DIR_FUNCTIONS


################################################################################
#########################         RandomForest         #########################
Feature_selection=T

### Train the hyperparameters ###
Rscript --vanilla $ML_method_rf $infile $Variable $PrevCutoff $Feature_selection $ncores $DIR_FUNCTIONS

### Stats best model ###
Dir_Results="./RandomForest_FS_TRUE/"
Rscript --vanilla $Plot_AUC_importance $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS

### Train the best model using the best hyperparameters ###
Rscript --vanilla $TrainModel_Complete_dataset $infile $Variable $Dir_Results $DIR_FUNCTIONS

################################################################################
#########################            XGBoost            ########################
#echo "XGBoost"
#Feature_selection=F
### The optimal one is the xgbGrid_mid
#GridMatrix="xgbGrid_mid"  #  GridMatrix <- "xgbGrid_large",GridMatrix <- "xgbGrid_mid",GridMatrix <- "xgbGrid_small", GridMatrix <- "xgbGrid_basic"

### Train the hyperparameters ###
#Rscript --vanilla $ML_method_XGboost $infile $Variable $PrevCutoff $Feature_selection $ncores $GridMatrix $DIR_FUNCTIONS

### Stats best model ###
#Dir_Results="./XGboost_FS_FALSE/"
#Rscript --vanilla $Plot_AUC_importance $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS

### Train the best model using the best hyperparameters ###
#Rscript --vanilla $TrainModel_Complete_dataset $infile $Variable $Dir_Results $DIR_FUNCTIONS




