set.seed(12345)
library(bnlearn)

argscomd = commandArgs(trailingOnly=TRUE) # argscomd <-"1"
IND <-  as.numeric(as.character(argscomd[1])) #1 # IND <- 1

path_func <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/functions"
source(paste0( path_func,"/Network_supplementary_functions.R" ))
source(paste0( path_func,"/bnlearn_supplementary_functions.R" ))


load("./boot_infiles/list_df.RData")

retBNlearns <- hc( list_df[[IND]], score = "custom", fun = my.bic.GLM, args = list() )
#retBNlearns <- hc( list_df[[IND]], score = "custom", fun = my.bic.GLM, args = list() , blacklist = Blist  )

dir.create("boot_outfiles")
saveRDS(file=paste0("./boot_outfiles/retBNlearns.",IND,".rds"))
