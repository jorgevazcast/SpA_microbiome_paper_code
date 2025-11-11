set.seed(12345)
library(bnlearn)
library(igraph)
library(bnlearn)

path_func <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code/functions"
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

