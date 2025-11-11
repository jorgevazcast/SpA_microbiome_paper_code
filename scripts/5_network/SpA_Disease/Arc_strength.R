set.seed(12345)
library(bnlearn)
library(parallel)

source("/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/functions/bnlearn_supplementary_functions.R")

df2bnlearn <- readRDS("./infile/df2bnlearn.rds")
nclusters <- 7

path_func <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/functions"

cl = makeCluster(nclusters)

# Exportar la variable 'path_func' al clúster
clusterExport(cl, varlist = c("path_func"))

# Cargar el archivo fuente dentro de cada worker
clusterEvalQ(cl, source("/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/functions/bnlearn_supplementary_functions.R"))

# Ejecutar el bootstrap
boot.strength.arcs = boot.strength(
  df2bnlearn,
  algorithm = "hc",
  algorithm.args = list(score = "custom", fun = my.bic.GLM),
  R = 500,
  cluster = cl
)

# Detener el clúster
stopCluster(cl)

# Guardar los resultados
dir.create("Arc_strength")
# Guardar los resultados
save(boot.strength.arcs, file = "./Arc_strength/boot.strength.arcs.RData")


