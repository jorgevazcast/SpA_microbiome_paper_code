set.seed(12345)
library(lmerTest)
working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

infile_metadata <- paste0(working_dir,"/Metadata/metadata_mice_exp.txt")
Metadata <- read.table(infile_metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(Metadata) <- gsub("-",".",Metadata$Barcode)
Metadata$Qfemto <- as.numeric(gsub(",",".",Metadata$Qfemto))
Metadata<-Metadata[complete.cases(Metadata$Qfemto),]
Metadata <- Metadata[Metadata$enterotype != "not_colonized",]

Metadata <- Metadata[!rownames(Metadata) %in% "Kinnex16S_Fwd_01..Kinnex16S_Rev_13",]
Metadata <- Metadata[!rownames(Metadata) %in% "X44_batch1",]

#### MEM #####
subdf <- subset(Metadata, Time == "w12" | Time == "w3" )
subdf$Sample_ID <- sapply(subdf$Sample,function(x){unlist(strsplit(x,"_"))[1]})
subdf$enterotype <- factor( subdf$enterotype, c("B2","B1") )

#formula_model <- paste(Y.var, "~" , paste(X.var,collapse = " + ")," + ", paste0("( 1 | ",Random_var, ")")   )
	
# rseModelSig <- lmer(Qfemto ~ Time + enterotype + sex +  Batch + (1|Sample_ID) , data= subdf )	
rseModelSig <- lmer(Qfemto ~ Time + enterotype + sex + (1|Batch) , data= subdf )	
	
sigModel <-  summary(rseModelSig) 
sigModel <- data.frame(sigModel$coefficients)
sigModel <-  data.frame(Var = rownames(sigModel) , sigModel )
sigModel$N <- nrow(subdf)	
sigModel$q.value <- sigModel$Pr...t..

write.table(sigModel, file = "MEM.tsv", sep = "\t", row.names = F)


