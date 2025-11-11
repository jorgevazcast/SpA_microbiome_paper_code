set.seed(12345)
library(bnlearn)
library(parallel)

working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

source(paste0(working_dir,"/functions/bnlearn_supplementary_functions.R"))

df2bnlearn <- readRDS("../infile/df2bnlearn.rds")
nclusters <- 7

path_func <- paste0(working_dir,"/functions")

cl = makeCluster(nclusters)

# Export the 'path_func' variable to all cluster nodes
clusterExport(cl, varlist = c("path_func"))

# Load the source file on each worker
clusterEvalQ(cl, source("~/github_shared_code_and_publications/SpA_microbiome_paper_code/functions/bnlearn_supplementary_functions.R"))

# Run the bootstrap procedure
boot.strength.arcs = boot.strength(
  df2bnlearn,
  algorithm = "hc",
  algorithm.args = list(score = "custom", fun = my.bic.GLM),
  R = 500,
  cluster = cl
)

# Stop the cluster
stopCluster(cl)

# Create directory for output
dir.create("Arc_strength")

# Save the results
save(boot.strength.arcs, file = "./Arc_strength/boot.strength.arcs.RData")

