set.seed(12345)
library(bnlearn)
library(igraph)

# Load custom functions
path_func <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/functions"
source(paste0(path_func, "/Network_supplementary_functions.R"))
source(paste0(path_func, "/bnlearn_supplementary_functions.R"))

# Load data
load("./boot_infiles/list_df.RData")

# Output directory
outdir <- "boot_outfiles"
dir.create(outdir, showWarnings = FALSE)

# Parallel processing
library(parallel)
#n_cores <- detectCores() - 1  # leave one core free
n_cores <- 10  # 

mclapply(1:length(list_df), function(IND) {
  cat("Running bootstrap ", IND, "\n")
  df <- list_df[[IND]]
  retBNlearns <- hc(df, score = "custom", fun = my.bic.GLM, args = list())
  saveRDS(retBNlearns, file = file.path(outdir, paste0("retBNlearns.", IND, ".rds")))
}, mc.cores = n_cores)




                       

    


