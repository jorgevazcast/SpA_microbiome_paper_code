# Bootstrap the dataset to generate resampled data for network stability analysis
Rscript --vanilla bootstrap_dataset.R  

# Submit the job for Bayesian network structure learning (Arc strength estimation) to the cluster
qsub Arc_strength.sh  

# Submit the job for Bayesian network bootstrapping using hill-climbing algorithm to the cluster
qsub Boot_bnlearn.sh  
