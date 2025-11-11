set.seed(12345)

df2bnlearn <- readRDS(file = "./infile/df2bnlearn.rds")

list_df <- list()
for(i in 1:500){
	subSamples <- sample(rownames(df2bnlearn), size=length(rownames(df2bnlearn)), replace = T )
	subDF <- df2bnlearn[match(subSamples,rownames(df2bnlearn)),]
	list_df[[i]] <- subDF
	rm(subSamples,subDF)
	
}

dir.create("boot_infiles")
save(file="./boot_infiles/list_df.RData",list_df)



