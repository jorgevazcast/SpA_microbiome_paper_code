# Rscript --vanilla Beta_diversity.R
set.seed(12345)
library(phyloseq)
library("ggplot2")
library(ggpubr)
library(gridExtra)

path_dir <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa"
source(paste0(path_dir, "/functions/beta_diver_functions_V9.R"))


PCoA_grapper <- function(Var = "", tempMetadata = data.frame(), bray.table.ADONIS=data.frame(),taxa_mat=matrix(),Distance="bray"){

	N_samples <- bray.table.ADONIS[Var,"N"]
	pval <-  round(  as.numeric(as.character(bray.table.ADONIS[Var,"p.value"])) ,digits=3)
	FDR <- round(  as.numeric(as.character(bray.table.ADONIS[Var,"BH.adj.p.value"])) ,digits=3)
	R2 <- round(  as.numeric(as.character(bray.table.ADONIS[Var,"R2"])) ,digits=3)	
	Fmodel <- round(  as.numeric(as.character(bray.table.ADONIS[Var,"Fmodel"])) ,digits=3)		

	
	PCoA_title <- paste0(Var,"\n"," p-val = ", pval,"; FDR = ", FDR,"; R2 = ",R2,"; F model = ", Fmodel,"; N = ", N_samples)

	tempMetadata$VarTemp <- tempMetadata[,Var]
	tempMetadata <- tempMetadata[!is.na(tempMetadata$VarTemp),]
	Levels <- sort(unique(tempMetadata$VarTemp))
	if(any(Levels == "HC")){Levels <- c("HC",Levels[!grepl("HC",Levels)])}
	
	if( length(sort(unique(tempMetadata$VarTemp))) == 2){colors<- c("blue","red");names(colors) <- Levels }
	if( length(sort(unique(tempMetadata$VarTemp))) == 3){colors<- c("blue","red","darkgreen");names(colors) <- Levels }
	if( length(sort(unique(tempMetadata$VarTemp))) == 4){colors<- c("blue","red","darkgreen","magenta");names(colors) <- Levels }
#	if( length(sort(unique(tempMetadata$VarTemp))) == 4){colors<- c("blue","#D95F02","#7570B3","#E7298A");names(colors) <- Levels }
#  
	taxa_mat <- taxa_mat[match( rownames(tempMetadata) , rownames(taxa_mat) ),]

	list_ordination <- vegan_PCoA_envfit(in.matrix= taxa_mat, Metadata2enfit=data.frame(tempMetadata), distance = Distance)
	PCoA<-list_ordination[["PCoA"]]; xlab<-list_ordination[["xlab"]]; ylab<-list_ordination[["ylab"]]; fit<-list_ordination[["fit"]]
	df_ord<-list_ordination[["df_ord"]]

	ret_list <- ggplot_envfit(df_ord, ord=PCoA,fit,Metadata=data.frame(tempMetadata),alpha_pval=0.1 )
		df_ord <- ret_list[["df_ord"]]
		df_arrows <- ret_list[["df_arrows"]]
		df_factors <- ret_list[["df_factors"]]
		
	if(Var == "Disease_activity"){
		colors <- c("#0000ee","#eead0e","#ee7621")
		df_ord$VarTemp <- factor(df_ord$VarTemp, levels=c("HC","SpALow", "SpAHigh"))
		PCoAplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab(xlab) + ylab(ylab)  + 
			stat_ellipse(aes(x = x, y = y,  group=VarTemp, fill=VarTemp ), 
			linetype = 2 ,type = "norm" , geom="polygon",level=0.8,alpha=0.1, show.legend=F) +	
			scale_fill_manual(values=colors ) + scale_color_manual(values=colors)+
			ggtitle(PCoA_title)  +  labs(color = Var)
			
	}else if(class(tempMetadata$VarTemp) == "character" | class(tempMetadata$VarTemp) == "factor" ){
	
		df_ord$VarTemp <- factor(as.character(df_ord$VarTemp), levels=Levels)
		PCoAplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab(xlab) + ylab(ylab)  + scale_color_manual(values=colors) + scale_fill_manual(values=colors) +
			stat_ellipse(aes(x = x, y = y,  group=VarTemp, fill=VarTemp ), 
			linetype = 2 ,type = "norm" , geom="polygon",level=0.8,alpha=0.1, show.legend=F) +	
			ggtitle(PCoA_title) +  labs(color = Var)
	}else{
		yourname_color_palette <- c("#74869c", "#6daddd", "#83adbb", "#b3d17c", "#ddaa7b","#ab5548")
		colors2use <- colorRampPalette(yourname_color_palette)(length(sort(unique(tempMetadata$VarTemp))))
		PCoAplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab(xlab) + ylab(ylab)  +	
			ggtitle(PCoA_title) +   scale_colour_gradientn(colors=colors2use) +  labs(color = Var) 	
	}
	return(PCoAplot)	
}


nmds_grapper <- function(Var = "", tempMetadata = data.frame(), bray.table.ADONIS=data.frame(),taxa_mat=matrix(),Distance="bray"){

	N_samples <- bray.table.ADONIS[Var,"N"]
	pval <-  round(  as.numeric(as.character(bray.table.ADONIS[Var,"p.value"])) ,digits=3)
	FDR <- round(  as.numeric(as.character(bray.table.ADONIS[Var,"BH.adj.p.value"])) ,digits=3)
	R2 <- round(  as.numeric(as.character(bray.table.ADONIS[Var,"R2"])) ,digits=3)	
	Fmodel <- round(  as.numeric(as.character(bray.table.ADONIS[Var,"Fmodel"])) ,digits=3)		

	
	NMDS_title <- paste0(Var,"\n"," p-val = ", pval,"; FDR = ", FDR,"; R2 = ",R2,"; F model = ", Fmodel,"; N = ", N_samples)

	tempMetadata$VarTemp <- tempMetadata[,Var]
	tempMetadata <- tempMetadata[!is.na(tempMetadata$VarTemp),]
	Levels <- sort(unique(tempMetadata$VarTemp))
	if(any(Levels == "HC")){Levels <- c("HC",Levels[!grepl("HC",Levels)])}
	
	if( length(sort(unique(tempMetadata$VarTemp))) == 2){colors<- c("blue","red");names(colors) <- Levels }
	if( length(sort(unique(tempMetadata$VarTemp))) == 3){colors<- c("blue","red","darkgreen");names(colors) <- Levels }
	if( length(sort(unique(tempMetadata$VarTemp))) == 4){colors<- c("blue","red","darkgreen","magenta");names(colors) <- Levels }
#	if( length(sort(unique(tempMetadata$VarTemp))) == 4){colors<- c("blue","#D95F02","#7570B3","#E7298A");names(colors) <- Levels }
#  
	taxa_mat <- taxa_mat[match( rownames(tempMetadata) , rownames(taxa_mat) ),]

	##### Perform a NMDS
	nmds_result <- metaMDS(taxa_mat, distance =Distance, k = 2, maxit = 999,  trymax = 500, wascores = TRUE)
	Caption <- paste0("Stress = ",round(nmds_result$stress,3))

#	goodness(nmds_result) # Produces a results of test statistics for goodness of fit for each point
#	stressplot(nmds_result) # Produces a Shepards diagram
	nmds_coordinates <- nmds_result$points

	# Create a data frame for ggplot2
	nmds_data <- as.data.frame(nmds_coordinates)
	colnames(nmds_data) <- c("x","y")
	nmds_data$Sample <- rownames(nmds_coordinates)
	nmds_data$VarTemp <- tempMetadata[match( rownames(nmds_data) ,  rownames(tempMetadata)),]$VarTemp
	df_ord <- nmds_data


	if(Var == "Disease_activity"){
		colors <- c("#0000ee","#eead0e","#ee7621")
		df_ord$VarTemp <- factor(df_ord$VarTemp, levels=c("HC","SpALow", "SpAHigh"))
		NMDSplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab("NMDS1") + ylab("NMDS2") + 
			stat_ellipse(aes(x = x, y = y,  group=VarTemp, fill=VarTemp ), 
			linetype = 2 ,type = "norm" , geom="polygon",level=0.8,alpha=0.1, show.legend=F) +	
			scale_fill_manual(values=colors ) + scale_color_manual(values=colors)+
			ggtitle(NMDS_title)  +  labs(color = Var, caption = Caption) 
			
	}else if(class(tempMetadata$VarTemp) == "character" | class(tempMetadata$VarTemp) == "factor" ){
	
		df_ord$VarTemp <- factor(as.character(df_ord$VarTemp), levels=Levels)
		NMDSplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab("NMDS1") + ylab("NMDS2") + scale_color_manual(values=colors) + scale_fill_manual(values=colors) +
			stat_ellipse(aes(x = x, y = y,  group=VarTemp, fill=VarTemp ), 
			linetype = 2 ,type = "norm" , geom="polygon",level=0.8,alpha=0.1, show.legend=F) +	
			ggtitle(NMDS_title) +  labs(color = Var, caption = Caption)
	}else{
		yourname_color_palette <- c("#74869c", "#6daddd", "#83adbb", "#b3d17c", "#ddaa7b","#ab5548")
		colors2use <- colorRampPalette(yourname_color_palette)(length(sort(unique(tempMetadata$VarTemp))))
		NMDSplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab("NMDS1") + ylab("NMDS2")  +	
			ggtitle(NMDS_title) +   scale_colour_gradientn(colors=colors2use) +  labs(color = Var, caption = Caption) 	
	}
	return(NMDSplot)	
}

	
