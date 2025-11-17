set.seed(12345)
library(phyloseq)
library(microbiome)

load_RData <- function(file) {
	# Load the data into a temporary environment
	temp_env <- new.env()
	load(file, envir = temp_env)
  
	# Get the first object (assuming there's only one)
	obj <- get(ls(temp_env)[1], envir = temp_env)
  
	# Return the renamed object
	return(obj)
}


read_input_phyloseq <- function(infile){

	if(grepl(".rds$",infile)){
		in.phylo <- readRDS(infile)
	}

	if( grepl(".RData$",infile) | grepl(".Rdata$",infile) ){
		in.phylo <- load_RData(file=infile)
	}
	return(in.phylo)

}


### The function create_phyloseq_object_fun can accepts matrix or a data.frame
create_phyloseq_object_fun <- function( data = data.frame() , tax = matrix(), sample_dat = data.frame()  ){
	cat("\n")	

	if (!is.matrix(tax)) {
		warning("The 'tax' object is not a matrix and will be converted to a matrix.")
		tax <- as.matrix(tax)
	}

	### Remove the X in the colnames of the samples, sometimes happens when the input dara colnames has 2.NAME as the first character of the name
	OTU = otu_table(data, taxa_are_rows = TRUE)
#	colnames(OTU)<-gsub("^X","",colnames(OTU))
	
	common_samples <- table(c(colnames(OTU) , rownames(sample_dat) ))
	common_samples <- names(common_samples[common_samples == 2])
	cat("Common samples = ", length(common_samples),"\n")
	cat("Exclusive samples from OTU table = ", length(setdiff( colnames(OTU) , common_samples  )),"\n")	
	cat("Exclusive samples from Metadata table = ", length(setdiff(  rownames(sample_dat) , common_samples  )),"\n")		
	
	common_taxa <- table(c(rownames(OTU) , rownames(tax) ))	
	common_taxa <- names(common_taxa[common_taxa == 2])
	cat("\nCommon taxa = ", length(common_taxa),"\n")
	cat("Exclusive taxa from OTU table = ", length(setdiff( rownames(OTU) , common_taxa  )),"\n")	
	cat("Exclusive taxa from Tax table = ", length(setdiff(  rownames(tax) , common_taxa  )),"\n")		
	
			
			
	OTU =  OTU[match(common_taxa,  rownames(OTU)),]
	OTU = OTU[order(rowSums(OTU),decreasing=T),]
	common_taxa <- rownames(OTU)
	TAX = tax_table(tax[match(common_taxa,  rownames(tax)),])
	
	SD = sample_data(sample_dat[match(common_samples,  rownames(sample_dat)),])
	
	OTU = otu_table( OTU[ , match(common_samples,  colnames(OTU))] , taxa_are_rows = T )	
	
	cat("\nN ", length(common_samples)," Congruent samples = ", all( sort(colnames(OTU)) == sort(rownames(SD))) , "\n")
	cat( "N ", length(common_taxa)," Congruent taxa = ", all( sort(rownames(OTU)) == sort(rownames(TAX))) , "\n")
	
	cat("\n")
	if( all( sort(colnames(OTU)) == sort(rownames(SD))) == FALSE ){
		print( head(  cbind(  sort(colnames(OTU)) ,  sort(rownames(SD)) )  ) )
		rownames(SD) <- gsub("[^[:alnum:] ]", ".", rownames(SD))
	}

	physeq = phyloseq(OTU, TAX,SD)
	return(physeq)
}

filter_and_match_data_phyloseq <- function(physeq.GMPR,physeq.rar,physeq.all, physeq.clr, inMetadata, Var2Use, prev_cutoff ){

	inMetadata$Var2USE <- inMetadata[[Var2Use]]
	inMetadata$SubjectID <- rownames(inMetadata)

	phyloForFilt <- create_phyloseq_object_fun( data = as.data.frame( as.matrix(physeq.rar@otu_table) ) , 
			tax = tax_table(physeq.rar), sample_dat = inMetadata  )
	inMetadata <- inMetadata[complete.cases(inMetadata$Var2USE),]
	phyloForFilt <- prune_samples(rownames(inMetadata) , phyloForFilt )
	Samples2use <- rownames(inMetadata)

	Taxa2keep <- c()
	for(i in unique(inMetadata$Var2USE)){
		
		# subphylo1 <- subset_samples(phyloForFilt, Var2USE== i)
		#subphylo <- subset_samples(phyloForFilt, eval(parse(text = paste0("Var2USE == '", i, "'"))) )
		samples_i <- rownames(inMetadata[inMetadata$Var2USE == i, ])
		subphylo <- prune_samples(samples_i, phyloForFilt)


		Prev <- prevalence(subphylo)
		Prev <- Prev[Prev >= prev_cutoff]
		Taxa2keep <- c( Taxa2keep , names(Prev))
		rm(Prev, subphylo,samples_i)
	}
	Taxa2keep <- unique(Taxa2keep)

	physeq.GMPR <-  prune_taxa(Taxa2keep, physeq.GMPR)
	physeq.rar <-  prune_taxa(Taxa2keep, physeq.rar)
	physeq.clr <-  prune_taxa(Taxa2keep, physeq.clr)
	physeq.all <-  prune_taxa(Taxa2keep, physeq.all)

	physeq.GMPR <-  prune_samples(Samples2use, physeq.GMPR)
	physeq.rar <-  prune_samples(Samples2use, physeq.rar)
	physeq.clr <-  prune_samples(Samples2use, physeq.clr)
	physeq.all <-  prune_samples(Samples2use, physeq.all)

	taxa_GMPR <-as.data.frame( as.matrix(physeq.GMPR@otu_table) )
	taxa_rar <-as.data.frame( as.matrix(physeq.rar@otu_table) )
	taxa_raw <-as.data.frame( as.matrix(physeq.all@otu_table) )
	taxa_clr <-as.data.frame( as.matrix(physeq.clr@otu_table) )

	taxa_GMPR <- taxa_GMPR[ , match( Samples2use ,  colnames(taxa_GMPR) ) ]
	taxa_raw <- taxa_raw[ , match( Samples2use ,  colnames(taxa_raw) ) ]
	taxa_rar <- taxa_rar[ , match( Samples2use ,  colnames(taxa_rar) ) ]
	taxa_clr <- taxa_clr[ , match( Samples2use ,  colnames(taxa_clr) ) ]

	rownames(taxa_GMPR) <- gsub("[^[:alnum:] ]", "_", rownames(taxa_GMPR)); rownames(taxa_GMPR) <- gsub(" ", "_", rownames(taxa_GMPR))
	rownames(taxa_raw) <- gsub("[^[:alnum:] ]", "_", rownames(taxa_raw)); rownames(taxa_raw) <- gsub(" ", "_", rownames(taxa_raw))
	rownames(taxa_rar) <- gsub("[^[:alnum:] ]", "_", rownames(taxa_rar)); rownames(taxa_rar) <- gsub(" ", "_", rownames(taxa_rar))
	Taxa2keep <- gsub("[^[:alnum:] ]", "_", Taxa2keep); Taxa2keep <- gsub(" ", "_", Taxa2keep)


#	taxa_GMPR <- taxa_GMPR[ match( Taxa2keep ,  rownames(taxa_GMPR) ) , ]
#	taxa_raw <- taxa_raw[ match( Taxa2keep ,  rownames(taxa_raw) ) , ]
#	taxa_rar <- taxa_rar[ match( Taxa2keep ,  rownames(taxa_rar) ) , ]
#	taxa_clr <- taxa_clr[ match( Taxa2keep ,  rownames(taxa_clr) ) , ]


	list_ret <- list(taxa_GMPR,taxa_raw,taxa_rar,taxa_clr, inMetadata)
	names(list_ret) <- c("taxa_GMPR","taxa_raw","taxa_rar","taxa_clr", "Metadata")
	return(list_ret)
}

filter_and_match_data_phyloseq_QMP <- function(physeq.rar, inMetadata, Var2Use, prev_cutoff ){
	inMetadata$Var2USE <- inMetadata[[Var2Use]]
	inMetadata$SubjectID <- rownames(inMetadata)
	
	if(is.null(physeq.rar@tax_table)){
		inTaxa <- as.matrix(data.frame( Taxa = taxa_names(physeq.rar), ID = taxa_names(physeq.rar)))
		rownames(inTaxa) <- taxa_names(physeq.rar)
		phyloForFilt <- create_phyloseq_object_fun( data = as.data.frame( as.matrix(physeq.rar@otu_table) ) , 
				tax = tax_table(inTaxa), sample_dat = inMetadata  )
					
	}else{
		phyloForFilt <- create_phyloseq_object_fun( data = as.data.frame( as.matrix(physeq.rar@otu_table) ) , 
				tax = tax_table(physeq.rar), sample_dat = inMetadata  )
	}

	inMetadata <- inMetadata[complete.cases(inMetadata$Var2USE),]
	phyloForFilt <- prune_samples(rownames(inMetadata) , phyloForFilt )
	Samples2use <- rownames(inMetadata)

	Taxa2keep <- c()
	for(i in unique(inMetadata$Var2USE)){
		
		# subphylo1 <- subset_samples(phyloForFilt, Var2USE== i)
		#subphylo <- subset_samples(phyloForFilt, eval(parse(text = paste0("Var2USE == '", i, "'"))) )
		samples_i <- rownames(inMetadata[inMetadata$Var2USE == i, ])
		subphylo <- prune_samples(samples_i, phyloForFilt)


		Prev <- prevalence(subphylo)
		Prev <- Prev[Prev >= prev_cutoff]
		Taxa2keep <- c( Taxa2keep , names(Prev))
		rm(Prev, subphylo,samples_i)
	}
	Taxa2keep <- unique(Taxa2keep)
	physeq.rar <-  prune_taxa(Taxa2keep, physeq.rar)	
	physeq.rar <-  prune_samples(Samples2use, physeq.rar)
	taxa_rar <-as.data.frame( as.matrix(physeq.rar@otu_table) )

	taxa_rar <- taxa_rar[ , match( Samples2use ,  colnames(taxa_rar) ) ]
	rownames(taxa_rar) <- gsub("[^[:alnum:] ]", "_", rownames(taxa_rar)); rownames(taxa_rar) <- gsub(" ", "_", rownames(taxa_rar))				
	#Taxa2keep <- gsub("[^[:alnum:] ]", "_", Taxa2keep); Taxa2keep <- gsub(" ", "_", Taxa2keep)
	list_ret <- list(taxa_rar, inMetadata)
	names(list_ret) <- c("taxa_rar", "Metadata")
	return(list_ret)	
}

filter_and_match_data_phyloseq_metabolome <- function(physeq.raw,physeq.norm, inMetadata, Var2Use, prev_cutoff ){
	inMetadata$Var2USE <- inMetadata[[Var2Use]]
	inMetadata$SubjectID <- rownames(inMetadata)
	
	if(is.null(physeq.raw@tax_table)){
		inTaxa <- as.matrix(data.frame( Taxa = taxa_names(physeq.raw), ID = taxa_names(physeq.raw)))
		rownames(inTaxa) <- taxa_names(physeq.raw)
		phyloForFilt <- create_phyloseq_object_fun( data = as.data.frame( as.matrix(physeq.raw@otu_table) ) , 
				tax = tax_table(inTaxa), sample_dat = inMetadata  )
					
	}else{
		phyloForFilt <- create_phyloseq_object_fun( data = as.data.frame( as.matrix(physeq.raw@otu_table) ) , 
				tax = tax_table(physeq.raw), sample_dat = inMetadata  )
	}

	inMetadata <- inMetadata[complete.cases(inMetadata$Var2USE),]
	phyloForFilt <- prune_samples(rownames(inMetadata) , phyloForFilt )
	Samples2use <- rownames(inMetadata)

	Taxa2keep <- c()
	for(i in unique(inMetadata$Var2USE)){
		
		# subphylo1 <- subset_samples(phyloForFilt, Var2USE== i)
		#subphylo <- subset_samples(phyloForFilt, eval(parse(text = paste0("Var2USE == '", i, "'"))) )
		samples_i <- rownames(inMetadata[inMetadata$Var2USE == i, ])
		subphylo <- prune_samples(samples_i, phyloForFilt)


		Prev <- prevalence(subphylo)
		Prev <- Prev[Prev >= prev_cutoff]
		Taxa2keep <- c( Taxa2keep , names(Prev))
		rm(Prev, subphylo,samples_i)
	}
	Taxa2keep <- unique(Taxa2keep)
	physeq.raw <-  prune_taxa(Taxa2keep, physeq.raw)	
	physeq.raw <-  prune_samples(Samples2use, physeq.raw)
	taxa_raw <-as.data.frame( as.matrix(physeq.raw@otu_table) )

	physeq.norm <-  prune_taxa(Taxa2keep, physeq.norm)	
	physeq.norm <-  prune_samples(Samples2use, physeq.norm)
	taxa_norm <-as.data.frame( as.matrix(physeq.norm@otu_table) )
	
	
	taxa_raw <- taxa_raw[ , match( Samples2use ,  colnames(taxa_raw) ) ]
	rownames(taxa_raw) <- gsub("[^[:alnum:] ]", "_", rownames(taxa_raw)); rownames(taxa_raw) <- gsub(" ", "_", rownames(taxa_raw))				
	#Taxa2keep <- gsub("[^[:alnum:] ]", "_", Taxa2keep); Taxa2keep <- gsub(" ", "_", Taxa2keep)

	taxa_norm <- taxa_norm[ , match( Samples2use ,  colnames(taxa_norm) ) ]
	rownames(taxa_norm) <- gsub("[^[:alnum:] ]", "_", rownames(taxa_norm)); rownames(taxa_norm) <- gsub(" ", "_", rownames(taxa_norm))				
	
	list_ret <- list(taxa_raw,taxa_norm, inMetadata)
	names(list_ret) <- c("taxa_raw", "taxa_norm", "Metadata")
	return(list_ret)	
}


contamination_taxa <- function(){
	#### Putative contaminants list
	# https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-014-0087-z
	contamination<-c("g_Afipia", "g_Aquabacterium", "g_Asticcacaulis", "g_Aurantimonas", "g_Beijerinckia", "g_Bosea", "g_Bradyrhizobium", "g_Brevundimonas", "g_Caulobacter", "g_Craurococcus", "g_Devosia", "g_Hoeflea", "g_Mesorhizobium", "g_Methylobacterium", "g_Novosphingobium", "g_Ochrobactrum", "g_Paracoccus", "g_Pedomicrobium", "g_Phyllobacterium", "g_Rhizobium", "g_Roseomonas", "g_Sphingobium", "g_Sphingomonas", "g_Sphingopyxis", "g_Acidovorax", "g_Azoarcus", "g_Azospira", "g_Burkholderiaf", "g_Comamonas", "g_Cupriavidus", "g_Curvibacter", "g_Delftia", "g_Duganella", "g_Herbaspirillum", "g_Janthinobacterium", "g_Kingella", "g_Leptothrix", "g_Limnobacter", "g_Massilia", "g_Methylophilus", "g_Methyloversatilis", "g_Oxalobacter", "g_Pelomonas", "g_Polaromonas", "g_Ralstonia", "g_Schlegelella", "g_Sulfuritalea", "g_Undibacterium", "g_Variovorax", "g_Acinetobacter", "g_Enhydrobacter", "g_Nevskia", "g_Pseudomonas", "g_Pseudoxanthomonas", "g_Psychrobacter", "g_Stenotrophomonas", "g_Xanthomonas", "g_Aeromicrobium", "g_Arthrobacter", "g_Beutenbergia", "g_Brevibacterium", "g_Corynebacterium", "g_Curtobacterium", "g_Dietzia", "g_Geodermatophilus", "g_Janibacter", "g_Kocuria", "g_Microbacterium", "g_Micrococcus", "g_Microlunatus", "g_Patulibacter", "g_Propionibacteriume", "g_Rhodococcus", "g_Tsukamurella", "g_Abiotrophia", "g_Bacillus", "g_Brevibacillus", "g_Brochothrix", "g_Facklamia", "g_Paenibacillus", "g_Streptococcus", "g_Chryseobacterium", "g_Dyadobacter", "g_Flavobacterium", "g_Hydrotalea", "g_Niastella", "g_Olivibacter", "g_Pedobacter", "g_Wautersiella", "g_Enterobacter", "g_Escherichia","g_Kineococcus","g_Ornithinimicrobium","g_Phenylobacterium","g_Escherichia/Shigella","g_Escherichia.Shigella","g_Tumebacillus","g_Escherichia_Shigella")

	### In case of having specific negative controls
	# CNsigcontamination<-read.table("/home/luna.kuleuven.be/u0120466/Postdoc_Raes/Projects/Prodigest/Analysis/Dysbiosis_and_deseases/Arthritis/SpA_paper_biopsies/Contamination/sv_cont_sig.tsv", header=T, row.names=1, dec=".", sep="\t")
	# more.cont <- rownames(subset(CNsigcontamination,Condition=="neg" & adj.pval< 0.1 ))
	# contamination<-unique(sort(c(contamination,more.cont)))
	return(contamination)
}
