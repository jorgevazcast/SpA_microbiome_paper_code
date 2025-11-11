set.seed(12345)
library(bnlearn)
library(igraph)
library(bnlearn)

#argscomd = commandArgs(trailingOnly=TRUE) # argscomd <-"1"
#IND <-  as.numeric(as.character(argscomd[1])) #1 # IND <- 1

path_func <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/functions"
source(paste0( path_func,"/Network_supplementary_functions.R" ))
source(paste0( path_func,"/bnlearn_supplementary_functions.R" ))

df2bnlearn <- readRDS(file="./infile/df2bnlearn.rds")



# retBNlearns <- hc( df2bnlearn, score = "custom", fun = my.bic.GLM, args = list() , blacklist = Blist  )
retBNlearns <- hc( df2bnlearn, score = "custom", fun = my.bic.GLM, args = list()  )

dir.create("outfile")
saveRDS(file="./outfile/retBNlearns.rds",retBNlearns)


g_bn <- as.igraph(retBNlearns)
pdf("./outfile/retBNlearns.pdf", width = 12, height =12)
	plot(g_bn, vertex.label.cex = 0.5, vertex.color = "lightblue",vertex.frame.color = "gray", 
		vertex.label.color = "black",edge.arrow.size = 0.5,layout = layout_with_fr,vertex.size = 8)
dev.off()

# mb(retBNlearns,"Linear.Tryptamine")
# mb(retBNlearns,"Boolean.enterotype.Bacteroides_2")
# mb(retBNlearns,"Boolean.Bacteroides_2")
# mb(retBNlearns,"Boolean.Disease.SpA")
# mb(retBNlearns,"Boolean.Disease_activity.SpAHigh")
# mb(retBNlearns,"Linear.Water")
# mb(retBNlearns,"Boolean.Disease_activity.HC")
# mb(retBNlearns,"QMP.Ruminococcus__gnavus__ref_mOTU_v25_01594_")

# mb(retBNlearns,"Linear.Fecal_calpro_values")
# mb(retBNlearns,"Boolean.Boolean.Disease_activity.SpAHigh")



                       

    


