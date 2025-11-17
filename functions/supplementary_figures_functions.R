library(rstatix)
library(ggmosaic)
library(ggstatsplot)
library(ggplot2)

working_dir <- "~/Postdoc_Raes/Projects/giant_cohort_spa" 
source(paste0(working_dir,"/functions/Extended_Statistical_Toolkit_functions.R"))
#path_estk <- "/home/luna.kuleuven.be/u0141268/github_projects/supplementary-statistical-functions"
#source(paste0(path_estk,"/Functions/Extended_Statistical_Toolkit_functions.R"))

color_varaibles <- function(Variable = ""){

	if(Variable == "Enterotype"){
		groupsOrder <- c("Bacteroides 1","Bacteroides 2","Prevotella","Ruminococcus")	
		colores <- c( "#f28118", "#ba1015", "#059571", "#4b59a1" )
		names(colores) <- groupsOrder
	}else if(Variable == "Bacteroides_2"){
		groupsOrder <- c("Bacteroides 2","Other")	
		colores <- c( "#ba1015", "gray" )
		names(colores) <- groupsOrder	
	}else if(Variable == "Bacteroides_2"){
		groupsOrder <- c("Bacteroides 2","Other")	
		colores <- c( "#ba1015", "gray" )
		names(colores) <- groupsOrder		
	}else if(Variable == "Disease_activity"){
		groupsOrder <- c("HC","SpALow","SpAHigh")
		colores <- c("#0000cd","#cdc673","#cd7054")
		names(colores) <- groupsOrder	
	}else if(Variable == "Diagnosis"){
		groupsOrder <- c("HC","Mixed","Axial","Peripheral")
		colores <- c("#0000ff","#cdc673","#7ccd7c","#cd3700")
		names(colores) <- groupsOrder	
	}else if(Variable == "Disease"){
		groupsOrder <- c("HC","SpA")
		colores <- c("blue","red")
		names(colores) <- groupsOrder			
	}else if(Variable == "Disease_activity"){
		groupsOrder <- c("HC","SpALow","SpAHigh")
		colores <- c("#0000ff","#cdc673","red")
		names(colores) <- groupsOrder
		
	}else if(Variable == "enterotype"){
		groupsOrder <- c("B1","B2")	
		colores <- c( "#f28118", "#ba1015" )
		names(colores) <- groupsOrder
	}else{
		colores <- NA
	}
		
	return(colores)
}



plot_stats_kw <- function(DiverMeta, colors2use, stats_df, Var, OrderFeat ){
	list_plots <- list()
	dunnTest_res <- data.frame()
	names_plots <- c()
	for(i in 1:nrow(stats_df)){

		temp_stas <- stats_df[i,]
		p_list <- cont_plot( subDF = DiverMeta,
			discrete = temp_stas$Discrete,
			continuous = temp_stas$Continuous,
			p.val = temp_stas$p.value,
			q.val = temp_stas$p.adjust,
			OrderFeat = names(colors2use),
			effectSize = temp_stas$Effsize_eta2,
			#Orderby = "median", # mean, median
			colors = colors2use,
			dunnTest = T)

		cat(temp_stas$Discrete,temp_stas$Continuous,"-- Kruskal–Wallis test:", 
			"p-value =",temp_stas$p.value, ", FDR =",temp_stas$p.adjust,
			", effect size (eta²) =",temp_stas$Effsize_eta2, ", N = ",temp_stas$Total_N,
			"\n")

		Plot <- p_list[[1]]
		Plot <- Plot + labs(fill = gsub("_"," ",Var))
		names_plots <- c(names_plots,temp_stas$Continuous)
	
		dunnTestRes <- p_list[[2]]
		colnames(dunnTestRes) <- gsub(".y.","Continuous",colnames(dunnTestRes))
		dunnTestRes$Continuous <- res_Diagnosis[i,]$Continuous
		dunnTest_res <- rbind(dunnTest_res, dunnTestRes)
		
		list_plots[[i]] <- Plot
		rm(temp_stas)
	}
	names(list_plots) <- names_plots
	dunnTest_res <- data.frame( Discrete = Var,  dunnTest_res)
	list_ret <- list(dunnTest_res, list_plots)
	names(list_ret) <- c("dunnTest_res","list_plots")
	return(list_ret)
}


plot_stats_wt <- function(DiverMeta, colors2use, stats_df, Var, OrderFeat ){
	list_plots <- list()
	dunnTest_res <- data.frame()
	names_plots <- c()
	for(i in 1:nrow(stats_df)){

		temp_stas <- stats_df[i,]
		p_list <- cont_plot( subDF = DiverMeta,
			discrete = temp_stas$Discrete,
			continuous = temp_stas$Continuous,
			p.val = temp_stas$p.value,
			q.val = temp_stas$p.adjust,
			OrderFeat = names(colors2use),
			effectSize = temp_stas$EffectSize_r,
			colors = colors2use,
			dunnTest = F)

		cat(temp_stas$Discrete,temp_stas$Continuous,"-- Wilcoxon signed-rank test:", 
			"p-value =",temp_stas$p.value, ", FDR =",temp_stas$p.adjust,
			", effect size (r) =",temp_stas$EffectSize_r, ", N = ",temp_stas$N,
			"\n")

		Plot <- p_list[[1]]
		Plot <- Plot + labs(fill = gsub("_"," ",Var))
		names_plots <- c(names_plots,temp_stas$Continuous)
		
		list_plots[[i]] <- Plot
		rm(temp_stas)
	}
	names(list_plots) <- names_plots
	return(list_plots)
}



