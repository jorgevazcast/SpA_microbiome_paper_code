set.seed(12345)
library(bnlearn)

path_func <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code/functions"

source(paste0( path_func,"/Network_supplementary_functions.R" ))
source(paste0( path_func,"/bnlearn_supplementary_functions.R" ))

argscomd = commandArgs(trailingOnly=TRUE) 
print(argscomd)

load("./boot_infiles/list_df.RData")

IndList <- as.numeric(as.character(argscomd[1])) #1
df2bnlearn <-list_df[[IndList]]
retBNlearns <- hc( df2bnlearn, score = "custom", fun = my.bic.GLM, args = list()  )
dir.create("out_networks")
save(file=paste0("./out_networks/retBNlearns.",IndList,".RData"),retBNlearns)
