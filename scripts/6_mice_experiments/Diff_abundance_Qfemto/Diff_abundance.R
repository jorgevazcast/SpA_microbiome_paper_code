set.seed(12345)
library("phyloseq")
library("ggplot2")
#library("microbiome")
library(lmerTest)
working_dir <- "~/github_shared_code_and_publications/SpA_microbiome_paper_code" 

source(paste0(working_dir,"/functions/Extended_Statistical_Toolkit_functions.R"))
source(paste0(working_dir,"/functions/supplementary_figures_functions.R"))

###############################################################################
############################       READ DATA       ############################
###############################################################################

#phyloseq_in <- readRDS(file=paste0(working_dir,"/1_infiles/mice_PacBio_data/phyloseq_20k_ee30/physeq_sv_Qfemto.rds"))
#Metadata <- data.frame(sample_data(phyloseq_in))
#Metadata<-Metadata[complete.cases(Metadata$Qfemto),]

infile_metadata <- paste0(working_dir,"/Metadata/metadata_mice_exp.txt")
Metadata <- read.table(infile_metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(Metadata) <- gsub("-",".",Metadata$Barcode)
Metadata$Qfemto <- as.numeric(gsub(",",".",Metadata$Qfemto))
Metadata<-Metadata[complete.cases(Metadata$Qfemto),]

Metadata <- Metadata[Metadata$enterotype != "not_colonized",]
Metadata <- Metadata[!rownames(Metadata) %in% "Kinnex16S_Fwd_01..Kinnex16S_Rev_13",]
Metadata <- Metadata[!rownames(Metadata) %in% "X44_batch1",]

binclass <- c("enterotype","sex")
multiclass <- c("cage","Extraction","Batch")

##############################
###### Test W3 binclass ######
subdf <- subset(Metadata, Time == "w3")
# fisher.test(table(subdf$enterotype,subdf$sex))
res_cond_bin <- data.frame()
for(i in binclass){
	tempRes <- cont_res_wt(subDF = subdf,discrete = i,continuous = "Qfemto")
	res_cond_bin <- rbind( res_cond_bin , tempRes )
}
res_cond_bin$p.adjust <-  p.adjust(res_cond_bin$p.value,method="BH")
write.table(res_cond_bin, file = "WT_res.w3.tsv", sep = "\t", row.names = F)

################################
###### Test W3 multiclass ######
res_cond_multi <- data.frame()
for(i in multiclass){
	tempRes <- cont_res_kw(subDF = subdf,discrete = i,continuous = "Qfemto")
	res_cond_multi <- rbind( res_cond_multi , tempRes[,c("Discrete", "Continuous","Kruskal.Wallis.chi.squared","Effsize_eta2","p.value","Total_N")] )

	
}
res_cond_multi$p.adjust <-  p.adjust(res_cond_multi$p.value,method="BH")
write.table(res_cond_multi, file = "KW_res.w3.tsv", sep = "\t", row.names = F)

##############################
###### Test W12 binclass ######
subdf <- subset(Metadata, Time == "w12")
# fisher.test(table(subdf$enterotype,subdf$sex))
res_cond <- data.frame()
for(i in binclass){
	tempRes <- cont_res_wt(subDF = subdf,discrete = i,continuous = "Qfemto")
	res_cond <- rbind( res_cond , tempRes )
}
res_cond$p.adjust <-  p.adjust(res_cond$p.value,method="BH")
write.table(res_cond, file = "WT_res.w12.tsv", sep = "\t", row.names = F)

################################
###### Test W12 multiclass ######
res_cond_multi <- data.frame()
for(i in multiclass){
	tempRes <- cont_res_kw(subDF = subdf,discrete = i,continuous = "Qfemto")
	res_cond_multi <- rbind( res_cond_multi , tempRes[,c("Discrete", "Continuous","Kruskal.Wallis.chi.squared","Effsize_eta2","p.value","Total_N")] )

	
}
res_cond_multi$p.adjust <-  p.adjust(res_cond_multi$p.value,method="BH")
write.table(res_cond_multi, file = "KW_res.w12.tsv", sep = "\t", row.names = F)

#######################
###### Time diff ######
subdf <- subset(Metadata, Time == "w12" | Time == "w3" )
subdf$Sample_ID <- sapply(subdf$Sample,function(x){unlist(strsplit(x,"_"))[1]})
subdf$enterotype <- factor( subdf$enterotype, c("B2","B1") )

M2U <- table(subdf$Sample_ID)
M2U <- sort(names(M2U[M2U == 2]))

w12 <- subset(subdf, Time == "w12"  )
w3 <- subset(subdf, Time == "w3"  )


w12 <- w12[match( M2U, w12$Sample_ID ),]
w3 <- w3[match( M2U, w3$Sample_ID ),]

indf <- rbind(w3,w12)

Time_diff <- cont_res_wt(subDF = indf,discrete = "Time",continuous = "Qfemto",Paired = T)
write.table(Time_diff, file = "WT_res.paired_w12_w3.tsv", sep = "\t", row.names = F)



