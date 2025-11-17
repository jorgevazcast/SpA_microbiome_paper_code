set.seed(12345)
library("phyloseq"); 
library("microbiome"); 

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
