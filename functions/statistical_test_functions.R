set.seed(12345)
library(ALDEx2)
library(DESeq2)
library(dunn.test)
library(GMPR)
library(coin)
library(ANCOMBC)
library(patchwork)
library(ACAT) # ACAT takes a series of p-values ​​(e.g., from different statistical tests applied to the same hypothesis or gene) and combines them into a single

#####################################################################
################      Wrapper for all the tests      ################
#####################################################################

RUN_tests <- function(Test, Meta, in_phylo, taxa.rar, taxa.log, taxa.raw, taxa.GMPR, taxa.clr, Var, QMPdata = F ){

 list_res <- vector(mode = "list", length = length(Test))
 names(list_res) <- Test
 
 if("ALDEx2" %in% Test){
	#########  ALDEx2   ##########
	res_ALDEx2 <- aldex_wrapper(in.meta = Meta, in.table = taxa.rar, Var2Use = Var )
	write.table(res_ALDEx2, "ALDEx2.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["ALDEx2"]] <- res_ALDEx2
 }
 if("DESeq2" %in% Test){
	#########  Deseq2   ##########
	res_Deseq2 <- deseq2_wrapper(in.meta = Meta, in.table = taxa.rar, Var2Use = Var, QMP = QMPdata )
	write.table(res_Deseq2, "DESeq2.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["DESeq2"]] <- res_Deseq2
 }
  if("GMPR_wt" %in% Test){
	#########  WT GMPR   ##########
	res_wt_GMPR <- diff_abund_wrapper(in.meta = Meta, in.table = taxa_GMPR, Var2Use = Var, method = "wt")
	write.table(res_wt_GMPR, "WT_GMPR.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["GMPR_wt"]] <- res_wt_GMPR		
 }
  if("RAR_wt" %in% Test){
	#########  WT RAR   ##########
	res_wt_rar <- diff_abund_wrapper(in.meta = Meta, in.table = taxa.rar, Var2Use = Var, method = "wt")
	write.table(res_wt_rar, "WT_rar.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["RAR_wt"]] <- res_wt_rar			
 }
  if("CLR_ttest" %in% Test){
	#########  Ttest CLR   ##########
	res_ttest_clr <- diff_abund_wrapper(in.meta = Meta, in.table = taxa.clr, Var2Use = Var, method = "t.test")
	write.table(res_ttest_clr, "Ttest_CLR.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["CLR_ttest"]] <- res_ttest_clr			
 }
  if("CLR_wt" %in% Test){
	#########  WT CLR   ##########
	res_wt_clr <- diff_abund_wrapper(in.meta = Meta, in.table = taxa.clr, Var2Use = Var, method = "wt")
	write.table(res_wt_clr, "WT_CLR.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["CLR_wt"]] <- res_wt_clr			
 } 

  if("log_ttest" %in% Test){
	#########  Ttest CLR   ##########
	res_ttest_clr <- diff_abund_wrapper(in.meta = Meta, in.table = taxa.log, Var2Use = Var, method = "t.test")
	write.table(res_ttest_clr, "Ttest_log.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["log_ttest"]] <- res_ttest_clr			
 }
 
 if("log_wt" %in% Test){
	#########  WT RAR   ##########
	res_wt_rar <- diff_abund_wrapper(in.meta = Meta, in.table = taxa.log, Var2Use = Var, method = "wt")
	write.table(res_wt_rar, "WT_log.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["log_wt"]] <- res_wt_rar			
 } 

 if("raw_wt" %in% Test){
	#########  WT RAR   ##########
	res_wt_rar <- diff_abund_wrapper(in.meta = Meta, in.table = taxa.raw, Var2Use = Var, method = "wt")
	write.table(res_wt_rar, "WT_raw.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["raw_wt"]] <- res_wt_rar			
 } 
  
 if("ANCOMBC" %in% Test){
	#########  ANCOMBC   ##########
	res_ancombc <- ancombc_wrapper(in.phylo = in_phylo, Var2Use = Var )
	write.table(res_ancombc, "ANCOMBC.results.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)
	list_res[["ANCOMBC"]] <- res_ancombc					
 }

 return(list_res)

}

#####################################################################
################           p-value matrix            ################
#####################################################################

join_pvalues_function <- function(Test_list, Ncat){


	if("ALDEx2" %in% names(Test_list)){
		if(Ncat == 2){
			aldex_names <- c("ALDEx2_ttest","ALDEx2_wt")
		}
	}

	###########################
	##### Common Features #####
	CommonFeatures <- c()
	for(i in names(Test_list)){
		print(i)
		Test_list[[i]]$Feature
		CommonFeatures <- c(CommonFeatures,Test_list[[i]]$Feature)
	}
	CommonFeatures <- table(CommonFeatures)
	CommonFeatures <- names(CommonFeatures[CommonFeatures == max(CommonFeatures)])

	###############################
	##### Join p-values table #####

	if("ALDEx2" %in% names(Test_list)){
		namescols <- names(Test_list)
		Ncol <- length(namescols) + 2
		p_values_table <- matrix(NA,length(CommonFeatures), Ncol ) 
		
		namescols <- namescols[!namescols %in% "ALDEx2"]
		namescols <- c("Feature",aldex_names,namescols)
		colnames(p_values_table) <- namescols
		p_values_table <- data.frame(p_values_table)
		p_values_table$Feature <- CommonFeatures

	}else{
		namescols <- names(Test_list)
		Ncol <- length(namescols) + 1
		p_values_table <- matrix(NA,length(CommonFeatures), Ncol )
		namescols <- c("Feature",namescols)
		colnames(p_values_table) <- namescols
		p_values_table <- data.frame(p_values_table)
		p_values_table$Feature <- CommonFeatures
		
	}


	#########################
	##### p-value table #####
	p_values_table_res <- fill_pval_matrix(in_list=Test_list,CF  = CommonFeatures,
				N = 2,in_tab = p_values_table,
				type="pval") # pval qval
	head(p_values_table_res)

	# # ACAT takes a series of p-values ​​(e.g., from different statistical tests applied to the same hypothesis or gene) and combines them into a single statistic using the Cauchy distribution.
	p_values_table_res$ACAT  <- apply(p_values_table_res[,2:ncol(p_values_table_res)], 1, function(pvec) ACAT(pvec))
	#p_values_table_res$ACAT  <- apply(p_values_table_res[,4:ncol(p_values_table_res)], 1, function(pvec) ACAT(pvec))

	p_values_table_res <- p_values_table_res[order( p_values_table_res$ACAT     ),]
	write.table(p_values_table_res, "p_values_table.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)

	#########################
	##### q-value table #####
	q_values_table <- p_values_table
	q_values_table_res <- fill_pval_matrix(in_list=Test_list,CF  = CommonFeatures,
				N = 2,in_tab = q_values_table,
				type="qval") # pval qval		
	
	ACAT.pval <- p_values_table_res$ACAT
	names(ACAT.pval) <- as.character(p_values_table_res$Feature)
	ACAT.pval <- ACAT.pval[match( CommonFeatures  , names(ACAT.pval) )]
	q_values_table_res$ACAT <- p.adjust(ACAT.pval,method="BH")

	q_values_table_res <- q_values_table_res[order( q_values_table_res$ACAT     ),]
	write.table(q_values_table_res, "q_values_table.tsv" ,quote=FALSE, sep= "\t",row.names = F, col.names = TRUE)

	list_ret <- list(p_values_table_res,q_values_table_res)
	names(list_ret) <- c("p_values_table_res","q_values_table_res")
	return(list_ret)

}


fill_pval_matrix <- function(in_list,CF,N,in_tab,type){  # type = "pval" # pval qval

	if(type == "pval"){
		if(N == 2){
			pval_lab <- "p.value" ; adelx2_ttest_pval_lab <- "we.ep"; adelx2_wt_pval_lab <- "wi.ep"			
		}
	}

	if(type == "qval"){
		if(N == 2){
			pval_lab <- "q.value" ; adelx2_ttest_pval_lab <- "we.eBH" ; adelx2_wt_pval_lab <- "wi.eBH"

		}				
	}
	for(i in names(in_list)){
		tempDF <- in_list[[i]]
		tempDF <- tempDF[match( CF  , tempDF$Feature ),]
		if(all(tempDF$Feature == in_tab$Feature)){
			 if( i  == "ALDEx2"){
			 	if(N == 2){
					in_tab[["ALDEx2_ttest"]] <- tempDF[[adelx2_ttest_pval_lab]]
					in_tab[["ALDEx2_wt"]] <- tempDF[[adelx2_wt_pval_lab]]
				}	 		
			 }else{
				in_tab[[i]] <- tempDF[[pval_lab]]
			 }
			rm(tempDF)
		}else{
			stop("Features do not match")
		}
	}

	return(in_tab)
}

##################################################################################
##################################### ALDEx2 #####################################
##################################################################################

#### Aldex wrapper funtion 
aldex_wrapper <- function(in.meta, in.table, Var2Use ){

	if( all( rownames(in.meta)  == colnames(in.table) ) ){
		print("Samples match")
	}else{
		in.table <- in.table[ , match( rownames(in.meta), colnames(in.table) ) ]
	}


	#conditions<-as.character(Metadata[colnames(in.sv.clr),Var2Use])
	conditions <- in.meta[[Var2Use]]
	names(conditions) <- rownames(in.meta)
	nr_var<-unique(sort(conditions))

	if( length(nr_var) ==  2){

		outfile<-paste("ALDEx2","tt",Var2Use,"tsv",sep=".")
		x <- aldex.clr(in.table, conditions, mc.samples=150, denom="clr", verbose=TRUE) # iqlr
	#	x.tt <- aldex.ttest(x, conditions, paired.test=FALSE)
		x.tt <- aldex.ttest(x, paired.test=FALSE)
		x.effect <- aldex.effect(x, include.sample.summary=FALSE, verbose=TRUE)
		ALDEx2.tt.x.all <- data.frame(x.tt,x.effect)
		ALDEx2.tt.x.all <- ALDEx2.tt.x.all[order(ALDEx2.tt.x.all$we.eBH),]
		#write.table(ALDEx2.tt.x.all, outfile,quote=FALSE, sep= "\t",row.names = T, col.names = TRUE) # Create the data table

		# Plot the pairwise ALDEx2 change
		ALDEx2plot(aldex2plot=ALDEx2.tt.x.all, adjusted_pvalue=0.1, num_cond = 2)

		ret_res <- ALDEx2.tt.x.all
		
	} else if( length(nr_var) >  2){
		#outfile<-paste("ALDEx2","glm",Var2Use,"tsv",sep=".")
		#outfile2<-paste("ALDEx2","pairwise_comparisons",Var2Use,"tsv",sep=".")
		x.all <- aldex(in.table, conditions, mc.samples=150, test="kw", effect=TRUE, include.sample.summary=TRUE, denom="all", verbose=FALSE)
		ALDEx2.glm.x.all <- data.frame(x.all)
		#write.table(ALDEx2.glm.x.all, outfile ,quote=FALSE, sep= "\t",row.names = T, col.names = TRUE) # Create the data table
		#x <- aldex.clr(in.table, conditions, mc.samples=150, denom="clr", verbose=TRUE)
		#x.glm <- aldex.glm(x, conditions)
		#ALDEx2.glm.x.all <- data.frame(x.glm)
		#write.table(ALDEx2.glm.x.all, outfile ,quote=FALSE, sep= "\t",row.names = T, col.names = TRUE) # Create the data table
		namesMetadata<-conditions
		ALDEx2.pairwsise<-sub_pairwise_aldex2(namesMetadata= namesMetadata, in.table = in.table)
		#write.table(ALDEx2.pairwsise, outfile2,quote=FALSE, sep= "\t",row.names = T, col.names = TRUE) # Create the data table

		# Plot the pairwise ALDEx2 change
		ALDEx2plot(aldex2plot=ALDEx2.pairwsise, adjusted_pvalue=0.1)
		ret_res <- ALDEx2.glm.x.all
				
	} else{

		stop("The variable contains 1 category")
	}
	ret_res <- data.frame(Feature = rownames(ret_res),ret_res)
	rownames(ret_res) <- NULL
	return(ret_res)
}



sub_pairwise_aldex2<- function(namesMetadata= "" , in.table = ""){

	### effect - median effect size:  diff.btw / max(dif.win) for all instances
	varMetadata <- unique(sort(namesMetadata))
	comb<-t(combn(varMetadata, 2))

	output<-matrix("", 0, 5  )
	colnames(output) <- c("comparison","Taxa","we.ep","we.eBH","effect")


	for(i in 1:dim(comb)[1]){
		print(i)
		cond1<-namesMetadata[ namesMetadata == comb[i,1] ]
		cond2<-namesMetadata[ namesMetadata == comb[i,2] ]
		### Subset ###
		sub.in.table<-in.table[,c(names(cond1),names(cond2))]
		Subconditions<-namesMetadata[colnames(sub.in.table)]

		ALDEx2.tt.x.all <- aldex(sub.in.table, Subconditions, mc.samples=150, denom="iqlr",  test="t", verbose=TRUE)
		comparison <-paste(comb[i,],collapse="-") 

	 
		ALDEx2.tt.x.all[,c("we.ep","we.eBH","effect")]
		tempmatrix<-cbind( rep(comparison,length(rownames(ALDEx2.tt.x.all))), 
				rownames(ALDEx2.tt.x.all) , ALDEx2.tt.x.all[,c("we.ep","we.eBH","effect")]  )
		colnames(tempmatrix) <- c("comparison","Taxa","we.ep","we.eBH","effect")

		output<-rbind(output,tempmatrix)

	}
	colnames(output) <- c("comparison","Taxa","we.ep","we.eBH","effect")
	return(output)
}


ALDEx2plot<-function(aldex2plot="", adjusted_pvalue=0.05, num_cond = ""){


	aldex2plot<-subset(aldex2plot, we.eBH < adjusted_pvalue)
	 if( dim(aldex2plot)[1] == 0 ){
		return(paste("Not significative features, P-value cut-off",adjusted_pvalue))
	} else if(num_cond == 2){


		aldex2plot$Effect <- ifelse(aldex2plot$effect < 0, "negative", "positive")  # above / below avg flag
		aldex2plot$Effect <- factor(aldex2plot$Effect,c("positive","negative"))

		aldex2plot<-aldex2plot[order(aldex2plot$effect,decreasing=F),]
		aldex2plot$Compare_categories <- rownames(aldex2plot)
		aldex2plot$Compare_categories <- factor(aldex2plot$Compare_categories,aldex2plot$Compare_categories)

		Num_features<-dim(aldex2plot)[1]
		pdfheight<-7
		# Diverging Barcharts
		deseqFig<-ggplot(aldex2plot, aes(x=Compare_categories, y=effect, label=Effect)) + theme_bw() +
			geom_bar(stat='identity', aes(fill=Effect), width=.5)  +
			scale_fill_manual(name="Effect change",  labels = c("Positive","Negative"),  values = c("positive"="#00ba38", "negative"="#f8766d")) + 
			labs(subtitle=paste("adjusted P value <",adjusted_pvalue),  title= "Pairwise ALDEx2 effect change") + 
			coord_flip()
		ggsave(paste0("ALDEx2_Effect_change",".pdf"),deseqFig,width = 11, height = pdfheight)
		
	
	} else{

		aldex2plot$Effect <- ifelse(aldex2plot$effect < 0, "negative", "positive")  # above / below avg flag
		aldex2plot$Effect <- factor(aldex2plot$Effect,c("positive","negative"))

		aldex2plot<-aldex2plot[order(aldex2plot$effect,decreasing=F),]
		aldex2plot$Compare_categories <- paste( aldex2plot$Taxa  ,  gsub("Condition-","",aldex2plot$comparison))
		aldex2plot$Compare_categories <- factor(aldex2plot$Compare_categories,aldex2plot$Compare_categories)

		Num_features<-dim(aldex2plot)[1]
		pdfheight<-7
		if(Num_features >= 30 &  Num_features < 35 ){pdfheight<-10}
		if(Num_features >= 50){pdfheight<-19}
		if(Num_features >= 100){pdfheight<-25}
		# Diverging Barcharts
		deseqFig<-ggplot(aldex2plot, aes(x=Compare_categories, y=effect, label=Effect)) + theme_bw() +
			geom_bar(stat='identity', aes(fill=Effect), width=.5)  +
			scale_fill_manual(name="Effect change",  labels = c("Positive","Negative"),  values = c("positive"="#00ba38", "negative"="#f8766d")) + 
			labs(subtitle=paste("adjusted P value <",adjusted_pvalue),  title= "Pairwise ALDEx2 effect change") + 
			coord_flip()
		ggsave(paste0("ALDEx2_Effect_change",".pdf"),deseqFig,width = 11, height = pdfheight)
      }
}


##################################################################################
##################################### DEseq2 #####################################
##################################################################################


deseq2_wrapper <- function(in.meta, in.table, Var2Use, QMP=F ){

	if( all( rownames(in.meta)  == colnames(in.table) ) ){
		print("Samples match")
	}else{
		in.table <- in.table[ , match( rownames(in.meta), colnames(in.table) ) ]
	}


	#conditions<-as.character(Metadata[colnames(in.sv.clr),Var2Use])
	conditions <- in.meta[[Var2Use]]
	names(conditions) <- rownames(in.meta)
	nr_var<-unique(sort(conditions))

	#### Create the metadata ####
	Metadata2test <- data.frame( Sample=names(conditions) , Condition=conditions)

	#### Estimate the size factor using GMPR  #####
	min_ct = 2 
	intersect_no = 4
	size.factor <- GMPR( (data.frame(t(in.table))) , min_ct = min_ct, intersect_no = intersect_no)
	size.factor <- size.factor[match( names(conditions), names(size.factor) )]

	#### Remove samples in case of NA data  #####
	if( any(is.na(size.factor)) == TRUE ){
		size.factor <- size.factor[!is.na(size.factor)]
		Metadata2test <- Metadata2test[ Metadata2test$Sample %in% names(size.factor),]
		in.table <- in.table[, match( names(size.factor)  , colnames(in.table) ) ]
		Metadata2test <- Metadata2test[ match( names(size.factor)  , Metadata2test$Sample ) ,  ]
	}

			
	if( length(nr_var) ==  2){

		#### Deseq analysis #####
		outfile<-paste("Deseq2",Var2Use,"tsv",sep=".")

		read.deseq2 <- DESeqDataSetFromMatrix(countData = in.table, colData = Metadata2test, design = ~Condition);
		cts <- counts(read.deseq2)
		geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
		dds.deseq2 <- estimateSizeFactors(read.deseq2, geoMeans=geoMeans)
		if(QMP==T){
			sizeFactors(dds.deseq2) <-  1
		}else{
			sizeFactors(dds.deseq2) <-  size.factor		
		}	
		#  Run deseq
		dds <- DESeq(dds.deseq2);
		res <- results(dds, test="Wald", contrast=c("Condition",nr_var[1],nr_var[2]), pAdjustMethod = "fdr", alpha = 0.05);
		res<- data.frame(res)
		res<- res[order(res$padj),]
		#write.table(res, outfile,quote=FALSE, sep= "\t",row.names = T, col.names = TRUE) # Create the data table

		# Plot the Log2 fold-change
		Deseq2plot(deseq2plot = res, adjusted_pvalue=0.1, num_cond = 2)
		
		ret_res <- res

	} else if( length(nr_var) >  2 ){

		#### Deseq analysis #####
		outfile<-paste("Deseq2","LRT",Var2Use,"tsv",sep=".")
		outfile2<-paste("Deseq2","pairwise_comparisons",Var2Use,"tsv",sep=".")

	# dds <- estimateSizeFactors(dds)

		in.table<-in.table[,names(conditions)]
		### LRT
		dds.LRT <- DESeqDataSetFromMatrix(countData = in.table, colData = Metadata2test, design = ~Condition)
		cts <- counts(dds.LRT)
		geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
		dds.LRT <- estimateSizeFactors(dds.LRT, geoMeans=geoMeans)
		if(QMP==T){
			sizeFactors(dds.LRT) <-  1
		}else{
			sizeFactors(dds.LRT) <-  size.factor		
		}	
		dds.LRT = DESeq(dds.LRT, test = "LRT", reduced = ~ 1)
		res.LRT<- data.frame(results(dds.LRT))
		res.LRT<- res.LRT[order(res.LRT$padj),]
		#write.table(res.LRT, outfile,quote=FALSE, sep= "\t",row.names = T, col.names = TRUE) # Create the data table

		### Pairwise Deseq2
		DeSeqInfoMatrix<-DESeqDataSetFromMatrix(countData = in.table, colData = Metadata2test, design = ~Condition)
		sizeFactors(DeSeqInfoMatrix) <-  size.factor
		res.LRT.pairwise<-sub_pairwise_deseq2( DeSeqInfoMatrix=DeSeqInfoMatrix, cond=conditions, condName="Condition" )
		res.LRT.pairwise<- res.LRT.pairwise[order(res.LRT.pairwise$padj),]
		write.table(res.LRT.pairwise, outfile2,quote=FALSE, sep= "\t",row.names = T, col.names = TRUE) # Create the data table

		# Plot the pairwise Log2 fold-change
		Deseq2plot(deseq2plot=res.LRT.pairwise, adjusted_pvalue=0.1)
	
		ret_res <- res.LRT

	} else{

		stop("The variable contains 1 category")
	}
	ret_res <- data.frame(Feature = rownames(ret_res),ret_res)
	rownames(ret_res) <- NULL
	colnames(ret_res) <- gsub("pvalue","p.value",colnames(ret_res))
	colnames(ret_res) <- gsub("padj","q.value",colnames(ret_res))
	return(ret_res)
}


##### Plot Deseq2 results
Deseq2plot<-function(deseq2plot="", adjusted_pvalue=0.05, num_cond = ""){


	deseq2plot<-subset(deseq2plot, padj < adjusted_pvalue)
	 if( dim(deseq2plot)[1] == 0 ){
		print(paste("Not significative features, P-value cut-off",adjusted_pvalue))
	}else if(num_cond == 2){
		deseq2plot$Log2FoldChange <- ifelse(deseq2plot$log2FoldChange < 0, "negative", "positive")  # above / below avg flag
		deseq2plot$Log2FoldChange <- factor(deseq2plot$Log2FoldChange,c("positive","negative"))

		deseq2plot$Compare_categories <- rownames(deseq2plot)
		deseq2plot<-deseq2plot[order(deseq2plot$log2FoldChange,decreasing=F),]
		deseq2plot$Compare_categories <- factor(deseq2plot$Compare_categories,deseq2plot$Compare_categories)
		Num_features<-dim(deseq2plot)[1]
			# Diverging Barcharts
		pdfheight<-7
		if(Num_features >= 30 &  Num_features < 35 ){pdfheight<-10}
		if(Num_features >= 50){pdfheight<-19}
		if(Num_features >= 100){pdfheight<-25}
		deseqFig<-ggplot(deseq2plot, aes(x=Compare_categories, y=log2FoldChange, label=Log2FoldChange)) + theme_bw() +
			geom_bar(stat='identity', aes(fill=Log2FoldChange), width=.5)  +
			scale_fill_manual(name="Fold Change",  labels = c("Positive","Negative"),  values = c("positive"="#00ba38", "negative"="#f8766d")) + 
			labs(subtitle=paste("adjusted P value <",adjusted_pvalue),  title= "Pairwise DESeq2  log2 fold change") + 
			coord_flip()
		ggsave(paste0("DESeq2_log2FoldChange",".pdf"),deseqFig,width = 11, height = pdfheight)

	}else{

	deseq2plot$Log2FoldChange <- ifelse(deseq2plot$log2FoldChange < 0, "negative", "positive")  # above / below avg flag
	deseq2plot$Log2FoldChange <- factor(deseq2plot$Log2FoldChange,c("positive","negative"))

	deseq2plot<-deseq2plot[order(deseq2plot$log2FoldChange,decreasing=F),]
	deseq2plot$Compare_categories <- paste( deseq2plot$Genus  ,  gsub("Condition-","",deseq2plot$comparison))
	deseq2plot$Compare_categories <- factor(deseq2plot$Compare_categories,deseq2plot$Compare_categories)

	Num_features<-dim(deseq2plot)[1]
	pdfheight<-7
	if(Num_features >= 30 &  Num_features < 35 ){pdfheight<-10}
	if(Num_features >= 50){pdfheight<-19}
	if(Num_features >= 100){pdfheight<-25}

	# Diverging Barcharts
	deseqFig<-ggplot(deseq2plot, aes(x=Compare_categories, y=log2FoldChange, label=Log2FoldChange)) + theme_bw() +
		geom_bar(stat='identity', aes(fill=Log2FoldChange), width=.5)  +
		scale_fill_manual(name="Fold Change",  labels = c("Positive","Negative"),  values = c("positive"="#00ba38", "negative"="#f8766d")) + 
		labs(subtitle=paste("adjusted P value <",adjusted_pvalue),  title= "Pairwise DESeq2  log2 fold change") + 
		coord_flip()
	ggsave(paste0("DESeq2_log2FoldChange",".pdf"),deseqFig,width = 11, height = pdfheight)
       }
}



#########################################################################################################
#####################################       Diff abund wrapper      #####################################
#########################################################################################################

diff_abund_wrapper <- function(in.meta, in.table, Var2Use, method = ""){ # wt t.test

	in.table <- data.frame(t(in.table))
	if( all( rownames(in.meta)  == rownames(in.table) ) ){
		print("Samples match")
	}else{
		in.table <- in.table[ match( rownames(in.meta), rownames(in.table) ) ,  ]
	}
	
	cont_variables <- colnames(in.table)
	in.table$Condition <- in.meta[[Var2Use]]

	ret_df <- data.frame()
	for(i in cont_variables){
		if(method == "wt"){
			temp_res <- cont_res_wt(subDF = in.table,discrete = "Condition",continuous = i)
		}
		if(method == "t.test"){
			temp_res <- cont_res_ttest(subDF = in.table,discrete = "Condition",continuous = i)		
		}		
		ret_df <- rbind( temp_res, ret_df )
		rm(temp_res)
	}
	
	ret_df$q.value <- p.adjust(ret_df$p.value, method = "BH")
	ret_df <- ret_df[order(ret_df$q.value),]
	ret_df$Discrete <- NULL
	colnames(ret_df) <- gsub("Continuous","Feature",colnames(ret_df))
	rownames(ret_df) <- NULL
	return(ret_df)
}



#########################################################################################################
#####################################   Wilcoxon signed-rank test   #####################################
#########################################################################################################


cont_res_wt <- function(subDF = data.frame(),discrete = "",continuous = ""){
	subDF <- subDF[,match(c(discrete,continuous), colnames(subDF))]
	colnames(subDF) <-c("discrete","continuous")
	subDF <- subDF[complete.cases(subDF),]
	subDF$discrete <- factor(subDF$discrete)
	
	WT <- wilcox.test(continuous ~ discrete, data = subDF)
	level1 = levels(subDF$discrete)[1]
	level2 = levels(subDF$discrete)[2]
	
	mean.level1 <-mean(subDF[subDF$discrete == level1,"continuous"])
	median.level1 <-median(subDF[subDF$discrete == level1,"continuous"])
	sd.level1 <-sd(subDF[subDF$discrete == level1,"continuous"])

	DFlevel1 <- subDF[subDF$discrete == level1,]
	DFlevel1 <- DFlevel1[complete.cases(DFlevel1),]
        N_level1 <- nrow(DFlevel1)
        
	mean.level2 <-mean(subDF[subDF$discrete == level2,"continuous"])
	median.level2 <-median(subDF[subDF$discrete == level2,"continuous"])
	sd.level2 <-sd(subDF[subDF$discrete == level2,"continuous"])  	

	DFlevel2 <- subDF[subDF$discrete == level2,]
	DFlevel2 <- DFlevel2[complete.cases(DFlevel2),]
        N_level2 <- nrow(DFlevel2)

	#####   Estimate the effect size	#####
	Z <- statistic(independence_test(continuous ~ discrete, data = subDF), type = "standardized")	
	N <- nrow(subDF)  # Total number of observations
	EffectSize_r <- c(abs(Z) / sqrt(N))
	names(EffectSize_r) <- "EffectSize_r"
	
	Dominant <- ifelse(median.level1 - median.level2 > 0, as.character(level1),as.character(level2))
	
	retDF<-data.frame(Discrete=discrete,Continuous=continuous,Level1=level1, N_level1,Level2=level2,N_level2,
		mean.level1,mean.level2,median.level1,median.level2,sd.level1,sd.level2,
		statistic=WT$statistic,EffectSize_r,p.value=WT$p.value, N=nrow(subDF),Dominant)
	return(retDF)
}



#########################################################################################################
#####################################        Student's t-test       #####################################
#########################################################################################################

cont_res_ttest <- function(subDF = data.frame(), discrete = "", continuous = "") {
	subDF <- subDF[,match(c(discrete,continuous), colnames(subDF))]
	colnames(subDF) <-c("discrete","continuous")
	subDF <- subDF[complete.cases(subDF),]
	subDF$discrete <- factor(subDF$discrete)

	if (nlevels(subDF$discrete) != 2) {
		stop("The variable 'discrete' must have exactly 2 levels for t-test.")
		}

	# t-test
	ttest <- t.test(continuous ~ discrete, data = subDF, var.equal = FALSE)

	level1 <- levels(subDF$discrete)[1]
	level2 <- levels(subDF$discrete)[2]

	# Summary statistics for level 1
	DFlevel1 <- subDF[subDF$discrete == level1, ]
	mean.level1 <- mean(DFlevel1$continuous)
	median.level1 <- median(DFlevel1$continuous)
	sd.level1 <- sd(DFlevel1$continuous)
	N_level1 <- nrow(DFlevel1)

	# Summary statistics for level 2
	DFlevel2 <- subDF[subDF$discrete == level2, ]
	mean.level2 <- mean(DFlevel2$continuous)
	median.level2 <- median(DFlevel2$continuous)
	sd.level2 <- sd(DFlevel2$continuous)
	N_level2 <- nrow(DFlevel2)

	# Cohen's d (tamaño del efecto)
	pooled_sd <- sqrt(((N_level1 - 1) * sd.level1^2 + (N_level2 - 1) * sd.level2^2) / (N_level1 + N_level2 - 2))
	cohen_d <- abs(mean.level1 - mean.level2) / pooled_sd
	names(cohen_d) <- "EffectSize_d"

	Dominant <- ifelse(median.level1 - median.level2 > 0, as.character(level1), as.character(level2))

	retDF <- data.frame( Discrete = discrete, Continuous = continuous, Level1 = level1, N_level1, Level2 = level2, N_level2,
		mean.level1, mean.level2, median.level1, median.level2, sd.level1, sd.level2,statistic = ttest$statistic,
		EffectSize_d = cohen_d, p.value = ttest$p.value,N = nrow(subDF),Dominant = Dominant
	)
	return(retDF)
}


#########################################################################################################
#####################################            ANCOMBC            #####################################
#########################################################################################################

ancombc_wrapper <- function(in.phylo, Var2Use ){

	out = ancombc(phyloseq = in.phylo, formula =Var2Use, p_adj_method = "BH", group = Var2Use, struc_zero = TRUE, neg_lb = FALSE,
                   tol = 1e-5, max_iter = 200, conserve = TRUE, alpha = 0.1, global = F) ## global T is for anova like results

	res <- combine_ancombc_results(out=out, filter_significant = FALSE, qval_cutoff = 0.1)
	rownames(res) <- NULL
	res <- res[order(res$q.value),]
	return(res)
}

combine_ancombc_results <- function(out, filter_significant = FALSE, qval_cutoff = 0.1) {

	# Extract result list
	res_list <- out$res

	# Convert individual result slots to data.frames
	dfs <- list(
		beta = as.data.frame(res_list$beta),
		se = as.data.frame(res_list$se),
		W = as.data.frame(res_list$W),
		pval = as.data.frame(res_list$p_val),
		qval = as.data.frame(res_list$q_val),
		diffabn = as.data.frame(res_list$diff_abn)
	)

	# Get reference rownames from one of the components
	taxa_names <- rownames(dfs$beta)

	# Check that all rownames are consistent across result data.frames
	consistent <- all(sapply(dfs, function(df) identical(rownames(df), taxa_names)))

	if (!consistent) {
		warning("Row names do not match across result tables. Please check manually.")

		# Show mismatches
		for (name in names(dfs)) {
			cat(paste0("Checking rownames in: ", name, "\n"))
			print(setdiff(rownames(dfs[[name]]), taxa_names))
		}

		return(NULL)
	}

	# Add suffix to column names to keep track of origin
	for (name in names(dfs)) {
		print(name)
		#colnames(dfs[[name]]) <- paste0(colnames(dfs[[name]]), "_", name)
		colnames(dfs[[name]]) <- name
	}

	# Combine all results
	res_df <- do.call(cbind, dfs)

	# Add Feature as column (from rownames)
	res_df$Feature <- rownames(res_df)
	res_df <- res_df[, c(ncol(res_df), 1:(ncol(res_df) - 1))]

	# Filter by significant results (qval < cutoff & diff_abn == TRUE), if requested
	if (filter_significant) {
		qval_cols <- grep("_qval$", colnames(res_df), value = TRUE)
		diffabn_cols <- grep("_diffabn$", colnames(res_df), value = TRUE)

		# Create logical mask per row if any qval < cutoff and diffabn == TRUE
		sig_mask <- apply(res_df[, qval_cols, drop = FALSE], 1, function(x) any(as.numeric(x) < qval_cutoff)) &
					apply(res_df[, diffabn_cols, drop = FALSE], 1, function(x) any(as.logical(x)))

		res_df <- res_df[sig_mask, , drop = FALSE]
	}
	colnames(res_df) <- gsub("pval","p.value",colnames(res_df))
	colnames(res_df) <- gsub("qval","q.value",colnames(res_df))	
	
	return(res_df)
}

#########################################################################################################
#####################################         Plot results          #####################################
#########################################################################################################


plot_significant_features <- function(in.meta, in.table, pval_table, Var2Use, cutoff = 0.1,colors=NULL, log_data = F ){

	Title_Line_plot <- "Effect Size (mean ± SD)"

	if( all( rownames(in.meta)  == colnames(in.table) ) ){
		print("Samples match")
	}else{
		in.table <- in.table[ , match( rownames(in.meta), colnames(in.table) ) ]
	}
	if(log_data == T){
		in.table <- log(in.table + 1)
		Title_Line_plot <- "Log Effect Size (mean ± SD)"
	}
	#conditions<-as.character(Metadata[colnames(in.sv.clr),Var2Use])
	conditions <- in.meta[[Var2Use]]
	names(conditions) <- rownames(in.meta)
	in.table <- data.frame(t(in.table))
	in.table <- in.table[ match(  names(conditions) , rownames(in.table) )  ,  ]
	
	in.table$Conditions <- conditions		

	in.table <- in.table[,colnames(in.table) %in%  c(as.character(pval_table$Feature),"Conditions")]

	summary_stats <- data.frame()
	for(i in unique(as.character(in.table$Conditions)) ){
		subdf <-  in.table[in.table$Conditions == i,]
		subdf$Conditions <- NULL
		Median <- apply(subdf, 2, median, na.rm = TRUE)
		Mean  <- apply(subdf, 2, mean, na.rm = TRUE)
		SD  <- apply(subdf, 2, sd, na.rm = TRUE)		
		Q1 <- apply(subdf, 2, function(x) quantile(x, probs = 0.25, na.rm = TRUE))
		Q3 <- apply(subdf, 2, function(x) quantile(x, probs = 0.75, na.rm = TRUE))

		summary_stats_temp <- data.frame( Condition = i, Feature = names(Median),Mean,SD,Median  = Median, Q1 = Q1, Q3 = Q3)
		summary_stats <- rbind(summary_stats,summary_stats_temp)
		rm(summary_stats_temp,subdf,Median,Mean,SD,Q1,Q3)

	}
	
	Feature_org <- names(sort(by( summary_stats$Mean, summary_stats$Feature ,  mean)))
	
	df_pvals_long <- reshape2::melt(pval_table)
	colnames(df_pvals_long) <- c("Feature","Test","q.value")
	df_pvals_long$Significant <- ifelse(df_pvals_long$q.value <= cutoff, "yes", "no")
	df_pvals_long$Feature <- factor(df_pvals_long$Feature, levels = Feature_org)
	
#	p1 <- ggplot(df_pvals_long, aes(x = Test, y = Feature, fill = Significant)) +
#		geom_tile(color = "white") +
#		scale_fill_manual(values = c("yes" = "black", "no" = "gray")) +
#		theme_bw() +
#		theme(
#			axis.text.x = element_text(angle = 45, hjust = 1),
#			legend.position = "bottom"
#		) +
#		ggtitle(paste0("Significance (p <= ",cutoff,")"))
	p1 <- ggplot(df_pvals_long, aes(x = Test, y = Feature, fill = Significant)) +
	  geom_tile(color = "black", linewidth = 0.5) +  # líneas negras entre tiles
	  geom_text(aes(label = round(q.value, 3)), color = "black", size = 3) +  # texto dentro de cada tile
	  scale_fill_manual(values = c("yes" = "gray", "no" = "white")) +
	  theme_bw() +
	  theme(
	    axis.text.x = element_text(angle = 45, hjust = 1),
	    legend.position = "bottom",
	    panel.grid = element_blank()  # elimina líneas de fondo
	  ) +
	  ggtitle(paste0("Significance (FDR <= ", cutoff, ")"))
	p1

	effect_df <- summary_stats
	effect_df$Feature <- factor(effect_df$Feature, levels = Feature_org)
		
	p2 <- ggplot(effect_df, aes(x = Mean, y = Feature, color = Condition)) +
	    geom_point(size = 3) +
	    geom_errorbarh(aes(xmin = Mean - SD, xmax = Mean + SD), height = 0.3) +

	    theme_bw() +
	    ggtitle(Title_Line_plot) +
	    theme(
	        axis.title.y = element_blank(),          # 
	        axis.text.y = element_blank(),           # 
	        axis.ticks.y = element_line(color = "black"),  # 
	        legend.position = "bottom"
	    )

	if( !is.null(colors) ){
		p2 <- p2 + scale_color_manual(values =colors[ levels(factor(effect_df$Condition))  ]) 
	}


	combined <- p1 + p2 + plot_layout(widths = c(0.7, 1.5))	  
	return(combined)	
}	



