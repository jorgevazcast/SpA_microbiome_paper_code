set.seed(12345)
library(phyloseq)
library(microbiome)
library(glmmTMB)
library(ggplot2)
library(dplyr)    # alternatively, this also loads %>%
#######################################################################################################
###########################################  Functions  ###############################################
#######################################################################################################

confound_variable <- function(taxa_sub = "", var_sub="", phy_df){

	### Subset the taxa and variable ####
	subdf <- phy_df[phy_df$OTU == taxa_sub,]
	subdf <- subdf[,c("Abundance",var_sub, "Water","BMI")]
	subdf <- subdf[complete.cases(subdf),]
	colnames(subdf) <- c("Abundance","VAR","Water","BMI")

	if(any(grepl("Other",subdf$VAR))){
	
		levelsVar<-unique(subdf$VAR)
		VarLevel <-levelsVar[!grepl("Other",levelsVar)]
		tempvar <- factor(as.character(subdf$VAR), levels=c("Other",var_sub) )
		subdf$VAR <- factor(subdf$VAR, levels=c("Other",var_sub) )
	}
#	print(range(subdf$Abundance))
	
	### Create the models ####	
	if(sum(subdf$Abundance == 0) / nrow(subdf) < 0.2){
		glm_nb.0 <- glmmTMB(formula = Abundance ~ Water + BMI , data = subdf,family =nbinom2 )
		glm_nb.1 <- glmmTMB(formula = Abundance ~ VAR + Water + BMI , data = subdf,family =nbinom2 )    
	}else{
		glm_nb.0 <- glmmTMB(formula = Abundance ~ Water + BMI , data = subdf,family =nbinom2, ziformula=~1 )
		glm_nb.1 <- glmmTMB(formula = Abundance ~ VAR + Water + BMI , data = subdf,family =nbinom2, ziformula=~1 )    

	}	
	
	# In case of the following message
	# Warning message:
	# In fitTMB(TMBStruc) :
  	#  Model convergence problem; extreme or very small eigenvalues detected. See vignette('troubleshooting')
  	# In some cases, extreme eigenvalues may be caused by having predictor variables that are on very different scales: try rescaling, and centering, continuous predictors in the model.
	na_return <- any(is.na(summary(glm_nb.1)$coefficients$cond[,4])) 
	if(na_return == TRUE){

		subdf$Water <- c(scale(subdf$Water))
		subdf$BMI <- c(scale(subdf$BMI))
		### Create the models ####	
		if(sum(subdf$Abundance == 0) / nrow(subdf) < 0.2){
			glm_nb.0 <- glmmTMB(formula = Abundance ~ Water + BMI, data = subdf,family =nbinom2 )
			glm_nb.1 <- glmmTMB(formula = Abundance ~ VAR + Water + BMI, data = subdf,family =nbinom2 )    
		}else{
			glm_nb.0 <- glmmTMB(formula = Abundance ~ Water + BMI, data = subdf,family =nbinom2, ziformula=~1 )
			glm_nb.1 <- glmmTMB(formula = Abundance ~ VAR + Water + BMI, data = subdf,family =nbinom2, ziformula=~1 )    

		}	
		
	
	}

	
	### Confounding variables ####	
	Anova_res <- data.frame(anova(glm_nb.0,glm_nb.1))
	
	### Return the data.frame  ####	
	Anova_res$Taxa <- taxa_sub
	Anova_res$VAR <- var_sub
	Anova_res$Model <- rownames(Anova_res)	
	
	#glm_nb.1.no_intercept <- update(glm_nb.1, ~ . -1) ### Show all the categories for the factor data

	glm1_coefficients <- data.frame(summary(glm_nb.1)$coefficients$cond)
	glm1_coefficients <- glm1_coefficients[grepl("VAR",rownames(glm1_coefficients)),]
	glm1_coefficients$Taxa <- taxa_sub
	glm1_coefficients$VAR <- var_sub
	glm1_coefficients$Model <- rownames(glm1_coefficients)	

	glm1_coefficients$var_cat <- paste0(glm1_coefficients$VAR,".",gsub("VAR","",glm1_coefficients$Model))

	
	list_res <- list(Anova_res,glm1_coefficients)	
	return(list_res)

	rm(Anova_res,glm_nb.1,glm_nb.0,glm_nb.1.no_intercept,subdf,glm1_coefficients)	
}


glmmTMB_grepper_model <- function(taxa_sub = "", var_sub="", phy_df){

	### Subset the taxa and variable ####
	subdf <- phy_df[phy_df$OTU == taxa_sub,]
	subdf <- subdf[,c("Abundance",var_sub, "Water","BMI")]
	subdf <- subdf[complete.cases(subdf),]
	colnames(subdf) <- c("Abundance","VAR","Water","BMI")

	if(any(grepl("Other",subdf$VAR))){
	
		levelsVar<-unique(subdf$VAR)
		VarLevel <-levelsVar[!grepl("Other",levelsVar)]
		subdf$VAR <- factor(subdf$VAR, levels=c("Other",var_sub) )
	}

	### Create the models ####	
	if(sum(subdf$Abundance == 0) / nrow(subdf) < 0.2){
		glm_nb.1 <- glmmTMB(formula = Abundance ~ VAR, data = subdf,family =nbinom2 )    
	}else{

		glm_nb.1 <- glmmTMB(formula = Abundance ~ VAR, data = subdf,family =nbinom2, ziformula=~1 )    

	}
	na_return <- any(is.na(summary(glm_nb.1)$coefficients$cond[,4])) 
	if(na_return == TRUE){
		subdf$Water <- c(scale(subdf$Water))
		subdf$BMI <- c(scale(subdf$BMI))
		### Create the models ####	
		if(sum(subdf$Abundance == 0) / nrow(subdf) < 0.2){
			glm_nb.1 <- glmmTMB(formula = Abundance ~ VAR, data = subdf,family =nbinom2 )    
		}else{
			glm_nb.1 <- glmmTMB(formula = Abundance ~ VAR, data = subdf,family =nbinom2, ziformula=~1 )    
		}
	}

	### Return the data.frame  ####	
	glm_coefficients <- data.frame(summary(glm_nb.1)$coefficients$cond)
	glm_coefficients$Taxa <- taxa_sub
	glm_coefficients$VAR <- var_sub
	glm_coefficients$Model <- rownames(glm_coefficients)	
		
	return(glm_coefficients)

	rm(glm_coefficients,glm_nb.1,subdf)	
}
corr_variable <- function(z=c(), meltPhyloseq=data.frame() ){
	x <- as.character(z[1])
	y <- as.character(z[2])
	xvec <- meltPhyloseq[ , colnames(meltPhyloseq) %in% x]
	yvec <- meltPhyloseq[ , colnames(meltPhyloseq) %in% y]
	bf<- cbind(xvec,yvec)
	bf<- bf[complete.cases(bf),]
	rcor <- cor.test(bf[,1], bf[,2],  method =  "spearman" )

	ret_df<- data.frame(statistic=rcor$statistic , estimate=rcor$estimate, p.value=rcor$p.value )
	ret_df$N <- length(bf[,1])
	return(ret_df)
}
rank_test_wrapper <- function(OTU_abundance, x_var){

	df2test<- data.frame( Abundance = OTU_abundance , Variable = x_var  )
	df2test<- df2test[complete.cases(df2test),]
	if(length(unique(df2test$Variable)) > 2){
		test_res <- kruskal.test(Abundance ~ Variable, data = df2test)
		rank_test_ret_df <- data.frame(statistic=test_res$statistic, p.value=test_res$p.value, method=test_res$method)
		rank_test_ret_df$N <- dim(df2test)[1]
	} else if( length(unique(df2test$Variable)) == 2 ){
		test_res <- wilcox.test(Abundance ~ Variable, data = df2test)
		rank_test_ret_df <- data.frame(statistic=test_res$statistic, p.value=test_res$p.value, method=test_res$method)
		rank_test_ret_df$N <- dim(df2test)[1]
	}
	return(rank_test_ret_df)

}
pair_test_variables_wrapper <- function(df_in = data.frame()){

	OTU_abundance <- df_in[,grepl("Abundance",colnames(df_in))]
	varaibles_categorical <- colnames(df_in); varaibles_categorical <- varaibles_categorical[!grepl("Abundance",varaibles_categorical)]
	ranktest_res<- lapply(varaibles_categorical, function(x_name,variables_data_frame=df_in, OTUAbundance=OTU_abundance ){
			#print(x_name)
			xvar <- variables_data_frame[,x_name==colnames(variables_data_frame) ]
			ret<-rank_test_wrapper(OTU_abundance=OTUAbundance, x_var=xvar)
			return(data.frame(variable=x_name,ret))
		})
	DF_cor <- do.call(rbind.data.frame, ranktest_res)
	rownames(DF_cor) <- NULL
	return(DF_cor)
}

binary_matrix <- function(var,var.name, return_factor = F){
	nr_var <- unique(var)
	nr_var <- nr_var[!is.na(nr_var)]
	matrix_ret<- matrix(var,length(var),length( nr_var ))
	matrix_ret[,1] == var
	colnames(matrix_ret) <- nr_var
	for(vanme in nr_var){
		matrix_ret[,vanme] <- ifelse(matrix_ret[,vanme] == vanme, paste0(var.name,".",vanme), "Other"   )
	}
	colnames(matrix_ret) <- paste0(var.name,".",nr_var)
	
	if(return_factor == T){
		
		for(i in 1:ncol(matrix_ret) ){
			matrix_ret[,i] <- ifelse(matrix_ret[,i] == "Other", 0, 1   )
		}
		matrix_ret <- data.frame(matrix_ret,stringsAsFactors=T)	
	}
	
	return(matrix_ret)
}


model.GLM = function(Y, X, data_in) {
	data_in <- data_in[ complete.cases(data_in[,c(Y,X)]) , ]
	data_in$Y <- data_in[[Y]]
	data_in$X <- data_in[[X]]	
	
	# Choose scoring strategy based on the type of node
	if (grepl("QMP", Y)) {
		ret <- model.nb( data = data_in) 	# Use negative binomial model
		sumret <-summary(ret)
		sumret <- sumret$coefficients$cond
	} else if (grepl("Linear", Y)) {
		ret <- model.lm( data = data_in) 	# Use linear model
		sumret <- summary(ret)
		sumret <- sumret$coefficients		
	} else if (grepl("Boolean", Y)) {
		ret <- model.logistic( data = data_in) 	# Use logistic regression
		sumret <- summary(ret)
		sumret <- sumret$coefficients				
	}
	return(sumret)
}

# Negative binomial BIC calculation using glmmTMB
model.nb = function( data) {

	formula_model <- as.formula("Y ~ X")
	# Create model formula depending on presence of parent nodes

	# Try fitting the model, use zero-inflation if needed
	ret <- tryCatch({
		if (sum(data[["Y"]] == 0) / nrow(data) < 0.2) {
			fit <- glmmTMB(formula = formula_model, data = data, family = nbinom2)
		} else {
			fit <- glmmTMB(formula = formula_model, data = data, family = nbinom2, ziformula = ~1)
		}
	}, error = function(e) NA) 	# Return NA if model fails

	return(fit)
}

# Linear model BIC calculation
model.lm = function( data) {
	formula_model <- as.formula("Y ~ X")
	# Fit linear model and return BIC
	ret <- tryCatch({
		fit <- lm(formula_model, data = data)
	}, error = function(e) NA) 	# Handle model failure

	return(fit)
}

# Logistic regression BIC calculation
model.logistic = function( data) {
	formula_model <- as.formula("Y ~ X")
	# Ensure binary variable has at least two levels
	if (length(unique(data[["Y"]])) < 2) {
		return(NA) 	# Cannot model with only one level
	}

	# Fit logistic model and return BIC
	ret <- tryCatch({
		fit <- glm(formula_model, family = binomial(link = 'logit'), data = data)
	}, error = function(e) NA) 	# Catch model errors

	return(fit)
}



