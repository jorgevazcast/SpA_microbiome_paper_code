set.seed(12345)
library(bnlearn)
library(igraph)

# Load custom functions
path_func <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/functions"
source(paste0(path_func, "/Network_supplementary_functions.R"))
source(paste0(path_func, "/bnlearn_supplementary_functions.R"))


#####################################################################################################################################################
##########################################################         FUNCTIONS             ############################################################
#####################################################################################################################################################

subset_boot_MB <-function(target_node="",List_files_boot_BN=list()){
	Target_MB <- c()
	for(f in List_files_boot_BN){
		load(f)
		retBNlearns$arcs
		Target_MB <- c(Target_MB ,c(mb(retBNlearns, target_node),target_node))
		rm(retBNlearns)
	}
	Target_MB <- data.frame(table(Target_MB))
	Target_MB<-Target_MB[order(Target_MB$Freq,decreasing=T),]
	Target_MB$Freq_norm<-Target_MB$Freq/length(List_files_boot_BN)
	return(Target_MB)

}

MB_network_igraph_old <- function(target_node = "", Boot_cutoff = 0.5, file_names, BayesianNetowkr){
	
	### Subset the MB ####
	vector.mb <- c(mb(retBNlearns, target_node),target_node)
	vector.mb  <- gsub("Linear.","",vector.mb )
	vector.mb  <- gsub("QMP.","",vector.mb )
	vector.mb  <- gsub("Boolean.","",vector.mb )

	### Subset the significant memebers of the MB (at least in 50% of the MB)
	Members.MB.BOOT <- subset_boot_MB(target_node=target_node,List_files_boot_BN=file_names)
	Members.MB.BOOT <- Members.MB.BOOT[Members.MB.BOOT$Freq_norm > Boot_cutoff,]
	Members.MB.BOOT$Target_MB  <- as.character(Members.MB.BOOT$Target_MB)
	Members.MB.BOOT$Target_MB  <- gsub("Linear.","",Members.MB.BOOT$Target_MB )
	Members.MB.BOOT$Target_MB  <- gsub("QMP.","",Members.MB.BOOT$Target_MB )
	Members.MB.BOOT$Target_MB  <- gsub("Boolean.","",Members.MB.BOOT$Target_MB )
	
	### Subset the MB
	vector.mb <- vector.mb[vector.mb %in% Members.MB.BOOT$Target_MB]
	subgraph <- induced_subgraph(complete_graph, vector.mb)
	size_vec  <- 20 * Members.MB.BOOT[match(V(subgraph)$name , Members.MB.BOOT$Target_MB ),]$Freq_norm
	V(subgraph)$size <- size_vec
	V(subgraph)$label.cex <- 1
	PercentBoot <- round(Members.MB.BOOT[match(V(subgraph)$name , Members.MB.BOOT$Target_MB ),]$Freq_norm,digits=2)
	V(subgraph)$name <- paste0( V(subgraph)$name , "\n Boot ", PercentBoot )

	return(subgraph)
}

markov_blanket <- function(NETWORK, target_node ){
	resMB <-  c(mb(retBNlearns, target_node),target_node) 
	return(resMB)
}		

MB_network_igraph <- function(target_node, ElemntsMB, NBootelements_MB, complete_graph, node_const_adapt = 30 ){
	subgraph <- induced_subgraph(complete_graph, ElemntsMB)
	Vsubname <- V(subgraph)$name
	BootNames <- paste( Vsubname , round(NBootelements_MB[Vsubname] / max(NBootelements_MB[Vsubname]) , digits =2 ) , sep = "\n") 
	V(subgraph)$name <- BootNames
	size_vec <-  (NBootelements_MB[Vsubname] / max(NBootelements_MB[Vsubname])) * node_const_adapt
	V(subgraph)$size <- size_vec
	plot(subgraph)

	return(subgraph)
}

		
#####################################################################################################################################################
##########################################################         SCRIPT             ###############################################################
#####################################################################################################################################################

############################################################
################      Read the infiles      ################

df2bnlearn <- readRDS("./infile/df2bnlearn.rds")

#####################################################################
################      Read the complete network      ################

retBNlearns <- readRDS("./outfile/retBNlearns.rds")

#############################################################################
################      Assing the coefficient for colors      ################
Arcs_cn <- data.frame(retBNlearns$arcs)
Estimates <- c()
for(i in 1:nrow(Arcs_cn)){
	print(i)
	subN <- Arcs_cn[i,]
	Estimates[i] <- model.GLM(Y=subN$from, X=subN$to, data_in=df2bnlearn)[2,"Estimate"] 
}
Arcs_cn$Estimates <- Estimates
Arcs_cn$color <- ifelse( Arcs_cn$Estimates > 0 , "blue" , "red"  )

##################################################################################
################      Assing the strenght and directionality      ################
load("./Arc_strength/boot.strength.arcs.RData")
StrDir_df <- data.frame()
for(i in 1:nrow(Arcs_cn)){
	subN <- Arcs_cn[i,]
	StrDir <- subset( boot.strength.arcs, from == subN$from & to ==  subN$to )
	StrDir_df <- rbind(StrDir_df,StrDir)	
}

Arcs_cn$strength <- StrDir_df$strength
Arcs_cn$direction <- StrDir_df$direction




#####################################################
###### Node type and color the phyloseq object ######
path_wd <- "/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/giant_cohort_spa/1_infiles"
##  physeq.qmp.motus
load(paste0(path_wd,"/QMP_mOTUS/physeq.qmp.motus.RData"))
taxa_names(physeq.qmp.motus) <- gsub("[^[:alnum:]]", "_", taxa_names(physeq.qmp.motus))
QMP_mOTUS <- gsub(" ", "_", taxa_names(physeq.qmp.motus))	
motus_names <- paste0("QMP.",QMP_mOTUS)

##  physeq.qmp.gmm
load(paste0(path_wd,"/QMP_GMM/physeq.qmp.gmm.RData"))
taxa_names(physeq.qmp.gmm) <- gsub("[^[:alnum:]]", "_", taxa_names(physeq.qmp.gmm))
QMP_GMM <- gsub(" ", "_", taxa_names(physeq.qmp.gmm))	
gmm_names <- paste0("QMP.",QMP_GMM)

##  physeq.metabolites.norm
load(paste0(path_wd,"/Metabolomics/physeq.metabolites.norm.RData"))
taxa_names(physeq.metabolites.norm) <- gsub("[^[:alnum:]]", "_", taxa_names(physeq.metabolites.norm))
metabolites <- gsub(" ", "_", taxa_names(physeq.metabolites.norm))	
metabolites_names <- paste0("Linear.",metabolites)

Metadata_var <- c("Boolean.Bacteroides_2","Boolean.Disease.SpA")

############################################################
################      Plot the network      ################
complete_graph <- graph_from_data_frame(Arcs_cn,directed = TRUE)
Vname <- V(complete_graph)$name

### Node colors
colors.nodes <- rep("gray",length(Vname))
names(colors.nodes) <- Vname
colors.nodes[names(colors.nodes) %in%  Metadata_var ] <- "red"
colors.nodes[names(colors.nodes) %in% gmm_names ] <- "#d67237"
colors.nodes[names(colors.nodes) %in% motus_names ] <- "darkgreen" # 
colors.nodes[names(colors.nodes) %in% metabolites_names ] <- "#5b1a18"
#colors.nodes[names(colors.nodes) %in% biomarkers_names ] <- "darkgreen" # 
colors.nodes.label <- colors.nodes
colors.nodes.plot  <- adjustcolor(colors.nodes, alpha = 0.6)
V(complete_graph)$color <- colors.nodes.plot
V(complete_graph)$label.color <- colors.nodes.label


### Frame color
frame.color  <- rep(NA,length(Vname))
V(complete_graph)$frame.color  <- frame.color 

#### Add the edge strenght and direction
E(complete_graph)$arrow.size <- Arcs_cn$strength 
E(complete_graph)$width <- Arcs_cn$strength * 5
E(complete_graph)$label = paste0("[S ",round(Arcs_cn$strength,digits=2),", D ", round(Arcs_cn$direction,digits=2),"]")
E(complete_graph)$label.cex = 1
E(complete_graph)$label.color = "black"
plot(complete_graph)

#### Plot the complete network
pdf("./outfile/complete_graph.pdf",width=15, height=12)
plot(complete_graph  , edge.curved=seq(-0.2, 0.2, length = ecount(complete_graph)) )
#plot(complete_graph    )
dev.off()

save(file="./outfile/complete_graph.RData", complete_graph)

############################################################
################      Markok blankets      ################

################      Count the the times that nodes are within the same Markov blanket    ################

# List all .rds files in the directory
files <- list.files("./boot_outfiles/", pattern = "\\.rds$", full.names = TRUE)

# Read them into a list
list_retBNlearns <- lapply(files, readRDS)

MB_freq_list <- vector("list", length(Vname))
names(MB_freq_list) <- Vname

for(target_node in names(MB_freq_list)){
	#target_node <- "QMP.Lachnospiraceae_species_incertae_sedis__meta_mOTU_v25_12240_"

	Target_MB <- c()
	for(i in 1:length(list_retBNlearns)){

		Target_MB <- c(Target_MB ,c(mb(list_retBNlearns[[i]], target_node) , target_node  ))
	}
	MB_freq_list[[target_node]] <- c(sort(table(Target_MB)))

}
MB_freq_list



################      Plot all the Markov blankets      ################

size_mb <- data.frame()
for(target_node in Vname ){
	MBtemp <- markov_blanket(NETWORK = retBNlearns, target_node =  target_node)
	MBlistfreq <-  MB_freq_list[[target_node]]
	
	mb_plot <- MB_network_igraph(target_node=target_node, 
		ElemntsMB=MBtemp, 
		NBootelements_MB=MBlistfreq, 
		complete_graph= complete_graph)
	
	print(target_node)
	print(length(MBtemp))
	size_mb <- rbind( size_mb ,data.frame( Node = target_node, N = length(MBtemp) ) )

	save(file=paste0("./outfile/MB.",target_node,".RData"), mb_plot)	
	pdf(file=paste0("./outfile/MB.",target_node,".pdf"),width=15, height=10)	
	 plot(mb_plot)
	 title(paste("MB", target_node))
	dev.off()
}

write.table(size_mb[order(size_mb$N,decreasing=T),], file = "./outfile/size_mb.tsv", col.names = TRUE, row.names = F,sep = "\t")


	
	
	
	
	
	


