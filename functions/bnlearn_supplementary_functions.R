set.seed(12345)
library(bnlearn)
library(glmmTMB)

### IMPORTANT ####
# https://www.bnlearn.com/documentation/man/score.html
# AIC and BIC are computed as logLik(x) - k * nparams(x), that is, the classic definition rescaled by -2. Therefore higher values are better, and for large sample sizes BIC converges to log(BDe).

#### RETURN THE NEGATIVE  BIC   #####################
# Main BIC scoring function used by bnlearn during structure learning
# Main BIC scoring function used by bnlearn during structure learning
my.bic.GLM = function(node, parents, data, args) {
	
	# Choose scoring strategy based on the type of node
	if (grepl("QMP", node)) {
		ret <- my.bic.nb(node, parents, data, args) 	# Use negative binomial model
	} else if (grepl("Linear", node)) {
		ret <- my.bic.lm(node, parents, data, args) 	# Use linear model
	} else if (grepl("Boolean", node)) {
		ret <- my.bic.logistic(node, parents, data, args) 	# Use logistic regression
	}

	# Handle failed model fits: return high penalty (Inf) if ret is NA or NULL
	if (is.null(ret) || is.na(ret)) {
		return(Inf)
	} else {
		return(-ret) 	# bnlearn maximizes score, so we negate the BIC
	}
}

# Negative binomial BIC calculation using glmmTMB
my.bic.nb = function(node, parents, data, args) {

	# Create model formula depending on presence of parent nodes
	if (length(parents) == 0) {
		model <- paste(node, "~ 1")
	} else {
		model <- paste(node, "~", paste(parents, collapse = "+"))
	}
	formula_model <- as.formula(model)

	# Try fitting the model, use zero-inflation if needed
	ret <- tryCatch({
		if (sum(data[[node]] == 0) / nrow(data) < 0.2) {
			fit <- glmmTMB(formula = formula_model, data = data, family = nbinom2)
		} else {
			fit <- glmmTMB(formula = formula_model, data = data, family = nbinom2, ziformula = ~1)
		}
		BIC(fit)
	}, error = function(e) NA) 	# Return NA if model fails

	return(ret)
}

# Linear model BIC calculation
my.bic.lm = function(node, parents, data, args) {

	# Create model formula
	if (length(parents) == 0) {
		model <- paste(node, "~ 1")
	} else {
		model <- paste(node, "~", paste(parents, collapse = "+"))
	}
	formula_model <- as.formula(model)

	# Fit linear model and return BIC
	ret <- tryCatch({
		fit <- lm(formula_model, data = data)
		BIC(fit)
	}, error = function(e) NA) 	# Handle model failure

	return(ret)
}

# Logistic regression BIC calculation
my.bic.logistic = function(node, parents, data, args) {

	# Create model formula
	if (length(parents) == 0) {
		model <- paste(node, "~ 1")
	} else {
		model <- paste(node, "~", paste(parents, collapse = "+"))
	}
	formula_model <- as.formula(model)

	# Ensure binary variable has at least two levels
	if (length(unique(data[[node]])) < 2) {
		return(NA) 	# Cannot model with only one level
	}

	# Fit logistic model and return BIC
	ret <- tryCatch({
		fit <- glm(formula_model, family = binomial(link = 'logit'), data = data)
		BIC(fit)
	}, error = function(e) NA) 	# Catch model errors

	return(ret)
}

