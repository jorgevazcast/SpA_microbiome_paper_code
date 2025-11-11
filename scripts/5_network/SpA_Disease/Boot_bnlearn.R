set.seed(12345)
library(bnlearn)

path_func <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/functions"

source(paste0( path_func,"/Network_supplementary_functions.R" ))
source(paste0( path_func,"/bnlearn_supplementary_functions.R" ))

argscomd = commandArgs(trailingOnly=TRUE) # argscomd <-"1"
print(argscomd)

#load("./boot_infiles/Blist.RData")
load("./boot_infiles/list_df.RData")
#Wlist <-  data.frame( from ="QMP.Ruminococcus__gnavus__ref_mOTU_v25_01594_" , to ="Linear.Tryptamine" )

IndList <- as.numeric(as.character(argscomd[1])) #1
df2bnlearn <-list_df[[IndList]]
#retBNlearns <- hc( df2bnlearn, score = "custom", fun = my.bic.GLM, args = list() , blacklist = Blist, whitelist= Wlist  )
retBNlearns <- hc( df2bnlearn, score = "custom", fun = my.bic.GLM, args = list()  )
dir.create("out_networks")
save(file=paste0("./out_networks/retBNlearns.",IndList,".RData"),retBNlearns)
#load("retBNlearns.RData")


