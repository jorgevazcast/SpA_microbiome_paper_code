set.seed(12345)
library(glmmTMB)    
library(glmnet)
library(lmerTest)
library(car)
####################################################################
##################       model functions          ##################
####################################################################

glmnet_wrapper <- function(Y = c(), X = matrix(), nfoldscv = 40,Alpha = 1, Nlambda=100){
	cvfit <- cv.glmnet(X,Y, nfolds = nfoldscv, nlambda = Nlambda, alpha = Alpha)
	Coeff <- coef(cvfit, s = "lambda.min")
	Coeff <- data.frame(as.matrix(Coeff))
	Coeff <- data.frame( Var = rownames(Coeff) , coeff.glmnet = Coeff$s1 )
	Coeff$Var <- gsub("[(]Intercept[)]","Intercept",Coeff$Var)
	return(Coeff)

}

glmnet_wrapper_table <- function(Group_variables_1 = c(), Group_variables_2 = c(), in_df = data.frame(), Nfoldscv = "LOOCV" ){ # Nfoldscv = 40 or "LOOCV"
	ret_glmnet <- data.frame()
	for(var in Group_variables_1){

		indata <- in_df[,colnames(in_df) %in% unique(c(Group_variables_2 , var))  ]
		indata <- indata[complete.cases(indata),]
		y = c(indata[,colnames(indata) %in% var])
		x = as.matrix(indata[,colnames(indata) %in% Group_variables_2])

		if(Nfoldscv == "LOOCV"){
			nfolds = length(y)
		}else{
			nfolds = Nfoldscv		
		}		
		df.ret.temp <- glmnet_wrapper(Y = y, X = x, nfoldscv = nfolds)
		df.ret.temp <- data.frame(Var1 = var,df.ret.temp, N = nrow(indata) )
		colnames(df.ret.temp) <- c("Var1","Var2","coeff.glmnet","N")
		ret_glmnet <- rbind(ret_glmnet,df.ret.temp)
		rm(indata,y,x,df.ret.temp)
	}
	ret_glmnet <- ret_glmnet[ret_glmnet$Var2 != "Intercept",]
	return(ret_glmnet)
}



nb_glm_function <- function(Y.var = "", X.var = c(""), in.df=data.frame()  ){
	formula_model <- paste(Y.var, "~" , paste(X.var,collapse = " + ") )
	#formula_model <- "Abundance ~ STUDY_GROUP + moisture + DEMOGRAPHIC_GENDER"
	formula_model <- as.formula(formula_model)

	#N <- nrow( in.df[complete.cases(in.df[,Y.var]), ]  )
	if( sum(in.df[,Y.var] == 0) / nrow(in.df) < 0.2 ){
		model_fit <- glmmTMB(formula = formula_model, data = in.df,family =nbinom2 )    
	}else{
		model_fit <- glmmTMB(formula = formula_model, data = in.df,family =nbinom2, ziformula=~1 )    
	}	
	
	Res_summ <- summary(model_fit)	
	Res_summ <- data.frame(Res_summ$coefficients$cond)
	categories = rownames(Res_summ) 
	var = rownames(Res_summ) 
	for(nam in X.var){ categories <- gsub(nam,"",categories) ; var[grepl(nam,var)] <- nam }
	categories[categories == ""] <- var[categories == ""]
	Res_summ <- data.frame( var, categories, Res_summ )
	rownames(Res_summ) <- NULL
	colnames(Res_summ)  <- gsub( "Pr...z.." , "p.value" , colnames(Res_summ) )
	colnames(Res_summ)  <- gsub( "Pr...t.." , "p.value" , colnames(Res_summ) )
	
	Res_anova <- data.frame(Anova(model_fit))	
	Res_anova <- data.frame( var = rownames(Res_anova) , Res_anova )	
	rownames(Res_anova) <- NULL

	list_ret <- list( Res_summ , Res_anova )
	names(list_ret) <- c("summary","anova")
	return(list_ret)
}

glm_function <- function(Y.var = "", X.var = c(""), in.df=data.frame()  ){
	formula_model <- paste(Y.var, "~" , paste(X.var,collapse = " + ") )
	#formula_model <- "Abundance ~ STUDY_GROUP + moisture + DEMOGRAPHIC_GENDER"
	formula_model <- as.formula(formula_model)
	
	model_fit <- glm(formula = formula_model, data = in.df )    
	
	Res_summ <- summary(model_fit)	
	Res_summ <- data.frame(Res_summ$coefficients)
	categories = rownames(Res_summ) 
	var = rownames(Res_summ) 
	for(nam in X.var){ categories <- gsub(nam,"",categories) ; var[grepl(nam,var)] <- nam }
	categories[categories == ""] <- var[categories == ""]
	Res_summ <- data.frame( var, categories, Res_summ )
	rownames(Res_summ) <- NULL
	colnames(Res_summ)  <- gsub( "Pr...z.." , "p.value" , colnames(Res_summ) )
	colnames(Res_summ)  <- gsub( "Pr...t.." , "p.value" , colnames(Res_summ) )


	Res_anova <- data.frame(Anova(model_fit))	
	Res_anova <- data.frame( var = rownames(Res_anova) , Res_anova )	
	rownames(Res_anova) <- NULL

	list_ret <- list( Res_summ , Res_anova )
	names(list_ret) <- c("summary","anova")
	return(list_ret)
}


### method = "lm"  method = "neg_binomial"
table_glm_test_function <- function(in.df=data.frame(), continuous_features = c(), Y.VAR = "", X.VAR = c() , 
			Feature = "OTU", Categorie2Explore = "", method = "lm" ){

	in.df$FEATURE = in.df[,Feature]
	ret_df <- data.frame()
	for(i in continuous_features){	
		sub_df <- subset(in.df, FEATURE == i)
		
		if(method == "neg_binomial"){	
			res_temp <- data.frame(Feature = i, nb_glm_function( Y.var = Y.VAR , X.var = X.VAR , in.df = sub_df )$summary )
		}
		if(method == "lm"){	
			res_temp <- data.frame(Feature = i, glm_function( Y.var = Y.VAR , X.var = X.VAR , in.df = sub_df )$summary )
		}			
		ret_df <- rbind(ret_df,res_temp )
		rm(sub_df,res_temp)
	}

	ret_df <- ret_df[ret_df$var == Categorie2Explore,]
	ret_df$q.val <- p.adjust(ret_df$p.value, method="BH") 
	
	if( any(colnames(ret_df) %in% "Pr...t..") ){ ret_df$q.val <- p.adjust(ret_df$Pr...t.., method="BH") }
	if( any(colnames(ret_df) %in% "Pr...z..") ){ ret_df$q.val <- p.adjust(ret_df$Pr...z.., method="BH")  }
		
	rownames(ret_df) <- NULL
	ret_df<-ret_df[order(ret_df$q.val),]
	return(ret_df)
}
##################################################################################
##################       Deconfounding model functions          ##################
##################################################################################

######################
#### Linear model ####
glm_cofound_function <- function(Y.var = "", NULL_model = c(""), Test_model = c(""), in.df = data.frame()) {

	# 1. Filter relevant columns and complete cases
	in.df <- in.df[, colnames(in.df) %in% unique(c(Y.var, NULL_model, Test_model))]
	in.df <- in.df[complete.cases(in.df), ]
	
	# Build model formulas
	formula_model_NULL <- as.formula(paste(Y.var, "~", paste(NULL_model, collapse = " + ")))
	formula_model_Test <- as.formula(paste(Y.var, "~", paste(Test_model, collapse = " + ")))

	# Fit models
	model_NULL <- glm(formula = formula_model_NULL, data = in.df)
	model_Test <- glm(formula = formula_model_Test, data = in.df)

	# Perform likelihood ratio test
	res_anova <- anova(model_NULL, model_Test, test = "Chisq")

	# Extract results (second row = comparison)
	  Res_df <- data.frame( Model = paste(Y.var, "~", paste(Test_model, collapse = " + ")), 
		  	Df = res_anova[2, "Df"], Deviance = res_anova[2, "Deviance"],Pr_Chi = res_anova[2, "Pr(>Chi)"])

	list_ret <- list(Res_df,model_Test)
	names(list_ret) <- c("Stats","Model")
	return(list_ret)
}

#################################
#### Negative binomial model ####
nb_cofound_function <- function(Y.var = "", NULL_model = c(""), Test_model = c(""), in.df = data.frame()) {

	# 1. Filter relevant columns and complete cases
	in.df <- in.df[, colnames(in.df) %in% unique(c(Y.var, NULL_model, Test_model))]
	in.df <- in.df[complete.cases(in.df), ]

	# Build model formulas
	formula_model_NULL <- as.formula(paste(Y.var, "~", paste(NULL_model, collapse = " + ")))
	formula_model_Test <- as.formula(paste(Y.var, "~", paste(Test_model, collapse = " + ")))
	# Fit models
	if( sum(in.df[,Y.var] == 0) / nrow(in.df) < 0.2 ){
		model_NULL <- glmmTMB(formula = formula_model_NULL, data = in.df,family =nbinom2 )
		model_Test <- glmmTMB(formula = formula_model_Test, data = in.df,family =nbinom2 )    		
	}else{
		model_NULL <- glmmTMB(formula = formula_model_NULL, data = in.df,family =nbinom2, ziformula=~1 )    
		model_Test <- glmmTMB(formula = formula_model_Test, data = in.df,family =nbinom2, ziformula=~1 )    
	}	
	# Perform likelihood ratio test
	res_anova <- anova(model_NULL, model_Test, test = "Chisq")
	# Extract results (second row = comparison)
	  Res_df <- data.frame( Model = paste(Y.var, "~", paste(Test_model, collapse = " + ")), 
		  	Df = res_anova[2, "Df"], Deviance = res_anova[2, "deviance"],Pr_Chi = res_anova[2, "Pr(>Chisq)"])

	list_ret <- list(Res_df,model_Test)
	names(list_ret) <- c("Stats","Model")
	return(list_ret)
}



#### Lemer ####
lmer_function_confound_variables <- function(Y.var = "", X.var = c(""), Random_var = "", in.df = data.frame()  ){

	#print(in.df)
	indf <<- in.df 
	
	formula_model <- paste(Y.var, "~" , paste(X.var,collapse = " + ")," + ", paste0("( 1 | ",Random_var, ")")   )
	formula_model <- as.formula(formula_model)

	rseModel <- lmer( formula_model, data= indf )

	cofound_res <- lmerTest::step(rseModel,direction = "both" )
	Fixed_res <- cofound_res$fixed
	
	if( any(Fixed_res$Eliminated == 0) ){
	
		X.var.sig <- rownames(Fixed_res[Fixed_res$Eliminated == 0,])
		formula_model_sig <- as.formula(paste(Y.var, "~" , paste(X.var.sig,collapse = " + ")," + ", paste0("( 1 | ",Random_var, ")") ))
		rseModelSig <- lmer( formula_model_sig, data= indf )		
		sigModel <- summary(rseModelSig)
		return(data.frame(sigModel$coefficients))
		#return(data.frame(Fixed_res[Fixed_res$Eliminated == 0,]))
	}else{
		return(NULL)
	}
}	


#################################
#### Step glm model ####
glm_step_function <- function(Y.var = "", Test_model = c(""), in.df = data.frame()) {

	# 1. Filter relevant columns and complete cases
	in.df <- in.df[, colnames(in.df) %in% c(Y.var, Test_model)]
	in.df <- in.df[complete.cases(in.df), ]

	# 2. Build the formula for the model
	formula_model <- as.formula(paste(Y.var, "~", paste(Test_model, collapse = " + ")))

	# 3. Fit the full model
	full_model <- glm(formula = formula_model, data = in.df)

	# 4. Select the best model based on AIC
	best_model <- step(full_model, direction = "both", trace = FALSE)

	# 5. Extract the summary and convert to data frame
	best_model_summary <- summary(best_model)
	best_model_coeff <- as.data.frame(best_model_summary$coefficients)

	SigVars <- as.character(best_model$formula)
	SigVars <- SigVars[SigVars %in% Test_model] 
	Res_summ <- best_model_coeff
	categories = rownames(Res_summ) 
	var = rownames(Res_summ) 
	for(nam in SigVars){ categories <- gsub(nam,"",categories) ; var[grepl(nam,var)] <- nam }
	categories[categories == ""] <- var[categories == ""]
	Res_summ <- data.frame( var, categories, Res_summ )
	rownames(Res_summ) <- NULL
	colnames(Res_summ)  <- gsub( "Pr...z.." , "p.value" , colnames(Res_summ) )
	colnames(Res_summ)  <- gsub( "Pr...t.." , "p.value" , colnames(Res_summ) )

	return(Res_summ)
}


