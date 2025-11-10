######################################################################################################
#########################         Varaible, cohort, and working dir         #########################
COHORT="Liver_cohort"
Variable="AST"
Normalization="Ret_abundance" ### Rarefaction  CLR GMPR

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
#HOME=/raeslab/scratch/jorvaz

#HOMEC=/raeslab/scratch/jorvaz
HOMEC=$HOME

DIR_FUNCTIONS=$HOMEC"/github_projects/mlpredictr/Functions" ## MODIFY THIS VARIABLE WHERE YOU HAVE THE REPO

#### Train the models and choose the ones with better hyperparameters
ML_method_rf=$HOMEC"/github_projects/mlpredictr/Scripts/RandomForest.R"
ML_method_glmnet=$HOMEC"/github_projects/mlpredictr/Scripts/glmnet.R"
ML_method_XGboost=$HOMEC"/github_projects/mlpredictr/Scripts/XGboost.R"

### Script for obtaining the stats for each model
Plot_RMSE_importance=$HOMEC"/github_projects/mlpredictr/Scripts/Plot_RMSE_importance.R"

### Train the complete dataset using the hyperparametesr of the best model
TrainModel_Complete_dataset=$HOMEC"/github_projects/mlpredictr/Scripts/TrainModel_Complete_dataset.R"

### Prevalence
PrevCutoff=0.2

### Ncores
ncores=10

####################################################################################
#########################         Infile                   #########################

infile="/home/luna.kuleuven.be/u0141268/Dropbox/Jiyeon_Jorge/research/research_liver_diseases/Predcition_ML/infiles/in_phylo_Genus_prop.rds"

################################################################################
#########################            GLMNET            #########################
#Feature_selection=F

### Train the hyperparameters ###
#Rscript --vanilla $ML_method_glmnet $infile $Variable $PrevCutoff $Feature_selection $ncores $DIR_FUNCTIONS

### Stats best model ###
#Dir_Results="./glmnet_FS_FALSE/"
#Rscript --vanilla $Plot_RMSE_importance $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS

### Train the best model using the best hyperparameters ###
#Rscript --vanilla $TrainModel_Complete_dataset $infile $Variable $Dir_Results $DIR_FUNCTIONS


################################################################################
#########################         RandomForest         #########################
Feature_selection=T

### Train the hyperparameters ###
Rscript --vanilla $ML_method_rf $infile $Variable $PrevCutoff $Feature_selection $ncores $DIR_FUNCTIONS

### Stats best model ###
Dir_Results="./RandomForest_FS_TRUE/"
Rscript --vanilla $Plot_RMSE_importance $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS

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
#Rscript --vanilla $Plot_RMSE_importance $infile $Variable $PrevCutoff $Dir_Results $DIR_FUNCTIONS

### Train the best model using the best hyperparameters ###
#Rscript --vanilla $TrainModel_Complete_dataset $infile $Variable $Dir_Results $DIR_FUNCTIONS




