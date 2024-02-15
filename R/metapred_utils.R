# Internal function for centering
# x numeric vector to be centered in data sets.
# center.in indices of data sets.
center <- function(x, center.in) {
  if (length(center.in) != length(x))
    stop("length(center.in) should match length(x).")
  for (trial in sort(unique(center.in))) {
    selection.id <- center.in == trial
    selection <- x[selection.id]
    x[selection.id] <- selection - mean(selection, na.rm = T)
  }
  x
}

# Center covariates within clusters
#
# Centers all except for y.name and cluster.name
#
# data data.frame
# y.name character, name of outcome variable
# cluster.name character, name of cluster variable.
centerCovs <- function(data, y.name, cluster.name) {
  to.center <- which((!(colnames(data) %in% cluster.name | colnames(data) %in% y.name) ) & sapply(data, is.numeric))
  cluster.vec <- data[ , cluster.name]
  
  for (col in to.center)
    data[ , col] <- center(data[ , col], center.in = cluster.vec)
  
  data
}

# These functions are used for making various names. Any change should be made here, such that
# these functions can also be used to retrieve objects from a list.
# st.u vector. unique indices or names of clusters.
# f numeric. fold index.
# type character. type of cv (just a chosen name)
# data.list list of data sets.
# data data.frame.
# covariate.columns indices of covariate columns.
# All return character.
# getFoldName <- function(st.u, f = NULL, type = NULL)
#   paste(getclName(st.u = st.u), getcvName(f = f, type = type), sep = if (length(f) > 0 || length(type) > 0) ". " else "")

# For l10, fixed and successive, the validation folds are unique. Bootstrap reuses them, and must add an iteration number.
getFoldName <- function(st.u, f = NULL, type = NULL) {
  if (isTRUE(type == "l1o") || isTRUE(type == "leaveOneOut") || isTRUE(type == "fixed") || isTRUE(type == "successive"))
    paste(getclName(st.u = st.u))
  else 
    paste(getclName(st.u = st.u), getcvName(f = f, type = type), sep = if (length(f) > 0 || length(type) > 0) ". " else "")
}

getclName <- function(st.u)
  paste(toString(st.u), sep = " ")

getcvName <- function(f, type = NULL)
  paste(type, f, sep = " ")

getStepName <- function(x)
  paste("s", x, sep = "")

getCoefs  <- function(fit, ...) {
  if (inherits(fit, "multinom"))
    return(coefMultinom(fit, ...))
  if (inherits(fit, c("glmerMod", "lmerMod")))
    return(lme4::fixef(fit, ...))
  else return(coef(fit))
}
# Needs some work. Should also return some coefficient names.
coefMultinom <- function(fit, ...)
  as.vector(t(coef(fit)))

# Perhaps unnecessary:
getVars   <- function(fit, ...) diag(as.matrix(vcov(fit)))
getCoVars <- function(fit, ...) as.matrix(vcov(fit))
getSE     <- function(fit, ...) sqrt(getVars(fit))

### The following functions are for generating the folst.u for the cross-validation in metapred
# st.u Numeric or character vector. Unique names of the strata / clusters.
# k Numeric. Differs per function
#   bootstrap: number of bootstraps
#   fixed: indices of data sets for validation.
#   successive (still a hidden function): Number of validation sets.
#   leaveOneOut = l1o: iecv: internal-external cross-validation of data sets/clusters.
leaveOneOut <- l1o <- function(st.u, ...) {
  st.u <- sort(st.u)
  if (length(st.u) < 2)
    stop("iecv not possible for fewer than 2 strata.")
  indexes <- seq_len(length(st.u))
  out <- list(dev = list(), dev.i = list(), val = as.list(st.u[indexes]), val.i = as.list(indexes))
  
  for (i in indexes) {
    out$dev[[i]]   <- st.u[-indexes[i]]
    out$dev.i[[i]] <- indexes[-i]
  }
  
  out
}

bs <- bootstrap <- function(st.u, k = NULL, ...) {
  st.u <- sort(st.u)
  if (is.null(k))
    k <- 200
  if (length(st.u) < 2)
    stop("Bootstrapping data sets is impossible for < 2 data sets.")
  dev <- dev.i <- val <- val.i <- list()
  
  i <- 1
  while (length(dev) < k) {
    indexes <- sample(1:length(st.u), replace = T)
    if (length(unique(indexes)) >= length(st.u))
      next
    dev[[i]] <- st.u[indexes]
    val[[i]] <- st.u[-indexes]
    dev.i[[i]] <- indexes
    val.i[[i]] <- seq_len(length(st.u))[-indexes]
    i <- i + 1
  }
  list(dev = dev, dev.i = dev.i , val = val, val.i = val.i)
}

fixed <- function(st.u, k = NULL, ...) {
  st.u <- sort(st.u)
  if (length(st.u) < 2)
    stop("Selecting a validation set is impossible for < 2 data sets.")
  if (is.null(k))
    k <- length(st.u)
  indexes <- seq_len(length(st.u))
  list(dev = list(st.u[-k]), dev.i = list(indexes[-k]), val = list(st.u[k]), val.i = list(indexes[k]))
}

ws <- within.sample <- function(st.u, ...) { # Necessary for testing recalibration functions.
  st.u <- sort(st.u)
  if (length(st.u) < 2)
    stop("Selecting a validation set is impossible for < 2 data sets.")
  indexes <- seq_len(length(st.u))
  list(dev = as.list(st.u[indexes]), dev.i = as.list(indexes), val = as.list(st.u[indexes]), val.i = as.list(indexes))
}

successive <- function(st.u, k = NULL, ...) {
  st.u <- sort(st.u)
  if (is.null(k)) k <- 1
  k <- as.integer(k)
  if (k < 1) stop("k must be >= 1")
  if (length(st.u) < (k + 1)) stop("Cross-validation requires k + 1 data sets, where k is the number of test sets.")
  out <- list(dev = list(), dev.i = list(), val = list(), val.i = list())
  
  for (i in seq_len(length(st.u) - k) ) {
    sel <- 1:i
    out$dev.i[[i]] <- sel
    out$dev[[i]] <- st.u[sel]
  }
  for (i in seq_len(length(out$dev)) ) {
    sel <- (i + 1):(i + k)
    out$val.i[[i]] <- sel
    out$val[[i]] <- st.u[sel]
  }
  out
}

# Remove observations (rows) that have at least one missing value.
# Necessary for metapred(), to ensure that the same observations are used
# df data.frame
# Returns: data.frame
remove.na.obs <- function(df) 
  df[apply(df, 1, function(df) sum(is.na(df))) == 0, ]

# Gets the predict method.
# fit Model fit object.
# two.stage logical. Is the model a two-stage model?
# predFUN Optional function, which is immediately returned
# ... For compatibility only.
getPredictMethod <- function(fit, two.stage = TRUE, predFUN = NULL, ...) {
  # A user written function may be supplied:
  if (!is.null(predFUN)) {
    if (is.function(predFUN)) {
      return(predFUN)
    } 
    return(get(as.character(predFUN), mode = "function"))
  }
  
  # If two-stage, the fit is used only to extract the link function.
  # If one-stage, fit's prediction method may be used.
  if (two.stage) {
    if (any(fit$stratified.fit[[1]]$stratum.class %in% c("logistf")) || inherits(fit, "logistf"))
      return(predictlogistf)
    if (any(fit$stratified.fit[[1]]$stratum.class %in% c("glm", "lm")) || inherits(fit, c("glm", "lm"))) 
      return(predictGLM)
    stop("No prediction method has been implemented for this model type yet for two-stage
         meta-analysis. You may supply one with the predFUN argument.")
  } else {
    if (any(fit$cv[[1]]$stratum.class %in% c("logistf")))
      return(predictlogistf)
    
    if (any(fit$cv[[1]]$stratum.class %in% c("glm", "lm")))
      return(predictGLM)
    
    if (any(fit$cv[[1]]$stratum.class %in% c("glmerMod", "lmerMod")))
      return(predictglmer)
  }
  
  # Return default predict function if everything else fails
  return(predict)
}

# Prediction function for two-stage metapred GLM objects
# object glm model fit object
# newdata newdata to predict for, of class "data.frame"
# b vector of coef. Overrides coef of object
# f formula used for selecting relevant variables from newdata. Overrides object
# ... For compatibility only.
# Returns vector of predicted values.
predictGLM <- function(object, newdata, b = NULL, f = NULL, type = "response", ...) {
  if (is.null(b)) b <- coef(object)
  if (is.null(f)) f <- formula(object)
  X <- model.matrix(f2rhsf(stats::as.formula(f)), data = newdata)
  
  lp <- X %*% b
  
  if (identical(type, "response")) {
    if (is.null(fam <- family(object)))
      return(lp)
    else
      return(fam$linkinv(lp))
  } else if (identical(type, "link"))
    return(lp)
}

predictglmer <- function(object, newdata, b = NULL, f = NULL, type = "response", ...)
  predictGLM(object = object, newdata = newdata, b = b, f = f, type = type, ...)

# Prediction function for logistf from the logisf package
# Args same as those of predictGLM()
predictlogistf <- function(object, newdata, b = NULL, f = NULL, type = "response", ...) {
  object$family <- binomial(link = "logit")
  predictGLM(object, newdata, b = b, f = f, type = type)
}

### NOTE: FOR SOME REASON UNBEKNOWN TO ME logistf() and requirenamespace(logistf) alter
# the functionalities of as.character(), thereby breaking the formula functions of
# this package.
### Note: cal.int does not work with this function. Use bin.cal.int instead.
# normal.int logical Should the intercept be recalibrated, such that Firth's correction
# is removed from it?
logistfirth <- function(formula = attr(data, "formula"), data = parent.frame(), pl = TRUE,
                        alpha = 0.05, control, plcontrol, firth = TRUE, init,
                        plconf = NULL, dataout = TRUE, ...) {
  if(is.null(normal.int <- list(...)$normal.int) ) normal.int <- FALSE
  if(is.null(fallback   <- list(...)$fallback)   ) fallback <- TRUE
  
  if (fallback && length(f2tl(formula)) < 1) 
    pl <- FALSE
  
  # Test if optional 'logistf' package is installed
  if (!requireNamespace("logistf", quietly = TRUE)) {
    stop("The package 'logistf' is currently not installed!")
  } 
  
  fit <- logistf::logistf(formula = formula, data = data, pl = pl,
                          alpha = alpha, control = control, plcontrol = plcontrol,
                          firth = firth, init = init,
                          plconf = plconf, dataout = dataout, ...)
  
  fit$family <- binomial()
  if (normal.int)
    return(fit)
  return(recalibrate(fit, newdata = data, f = ~ 1, estFUN = glm, family = binomial))
}

# Univariate Random Effects Meta-Analysis
# coefficients data.frame or matrix, containing coef
# variances data.frame or matrix, containing variances
# method Method for meta-analysis.
# vcov Ignored. For compatibility only, until a better solution is implemented.
# ... Optional arguments for rma().
#' @importFrom metafor rma
urma <- function(coefficients, variances, method = "DL", vcov = NULL, ...) {
  if (!(is.data.frame(coefficients) || is.matrix(coefficients)) || !(is.data.frame(variances) || is.matrix(variances)) )
    stop("coefficients and variances must both be a data.frame or matrix.")
  if (!identical(dim(coefficients), dim(variances)))
    stop("coefficients and variances must have the same dimensions.")
  
  meta.b <- meta.se <- meta.tau2 <- meta.ci.lb <- meta.ci.ub <- meta.pi.lb <- meta.pi.ub <-
    meta.se.tau2 <- rep(NA, ncol(coefficients))
  
  for (col in seq_len(ncol(coefficients))) {
    tryCatch(
      r <- metafor::rma(coefficients[ , col] , variances[ , col], method = method, ...),
      error = function(e) {
        stop(paste("Error in univariate rma of variable:", names(coefficients)[col], ",as follows:", e))
      }
    )
    
    meta.b[col]     <- r$beta
    meta.se[col]    <- r$se
    meta.tau2[col]  <- r$tau2
    meta.ci.lb[col] <- r$ci.lb
    meta.ci.ub[col] <- r$ci.ub
    
    # Warning because this is highly unexpected!
    # rma may change to fixed effects under unknown conditions (probably low sample size/ low number of clusters)
    # checks because otherwise an error is returned.
    if (!identical(r$method, method))
      warning(paste("metafor::rma switched automatically from ", method, " to ", r$method, " method.", sep = ""))
    
    if (!identical(r$method, "FE")) {
      cr <- metafor::predict.rma(r)
      meta.pi.lb[col] <- cr$cr.lb
      meta.pi.ub[col] <- cr$cr.ub
      meta.se.tau2[col] <- r$se.tau2
    }
  }
  
  meta.v <- meta.se^2
  names(meta.b) <- names(meta.v) <- names(meta.se) <- names(meta.tau2) <-  names(meta.ci.lb) <- 
    names(meta.ci.ub) <- names(meta.pi.lb) <-  names(meta.pi.ub) <- colnames(coefficients)
  
  list(coefficients = meta.b, variances = meta.v, se = meta.se, ci.lb = meta.ci.lb, ci.ub = meta.ci.ub,
       tau2 = meta.tau2, se.tau2 = meta.se.tau2, tau = sqrt(meta.tau2), 
       pi.lb = meta.pi.lb, pi.ub = meta.pi.ub, method = r$method)
}

# Multivariate Random Effects Meta-Analysis # Old function, see new function below.
# coefficients data.frame or matrix, containing coef
# vcov vcov as in vcov(stratified.fit)
# variances IGNORED.
# TBI: method Method for meta-analysis.
# TBI: ... Optional arguments for mvmeta().
# #' @importFrom mvmeta mvmeta
# mrma <- function(coefficients, vcov, variances, ...) {
#   # Test if optional 'mvmeta' package is installed
#   if (!requireNamespace("mvmeta", quietly = TRUE))
#     stop("The package 'mvmeta' is currently not installed.")
#   
#   # mvmeta assumes a differnt form for vcov when nvar = 1.
#   one_var <- identical(dim(coefficients)[2], 1L)
#   if (one_var)
#     vcov <- as.matrix(vcov)
#   
#   # fit   
#   ma.fit <- mvmeta::mvmeta(as.matrix(coefficients), S = vcov)
#   
#   # Rename, as mvmeta changes the names.
#   vcov <- ma.fit$vcov
#   colnames(vcov) <- row.names(vcov) <- colnames(coefficients)
#   
#   if (one_var)
#     coefficients <- ma.fit$coefficients[1]
#   else
#     coefficients <- ma.fit$coefficients[1,] # [1,] to coerce to vector with names.
#   
#   # Return
#   list(coefficients = coefficients, 
#        variances = diag(vcov),
#        se = sqrt(diag(vcov)),
#        vcov = vcov,
#        psi = ma.fit$Psi)
# }

# What is this?
blockMatrixDiagonal <- function(...){  
  matrixList <- list(...)
  if(is.list(matrixList[[1]])) matrixList <- matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
}

## vcov is a 3-dimensional array, where [,,i] indicates the within-study covariance of study i
mrma <- function(coefficients, vcov, variances,  ...) {
  meta.method <- "REML" # Method for meta-analysis
  
  # Test if we are dealing with multivariate data
  if (inherits(coefficients, "numeric")) {
    coefficients <- data.frame(y=coefficients)
  }
  if (ncol(coefficients)==1) {
    S = data.frame(S=vcov)
    return(urma(coefficients = coefficients, variances = S, method = meta.method))
  }
  
  # In rma.mv, we need to supply a stacked version of all coefficients, together with a grouping variable
  yi <- as.vector(t(coefficients))
  
  # Outcome indicator
  group <- rep(NA, length(yi))
  for (i in 1:ncol(coefficients)) {
    group[seq(i,length(yi), by=ncol(coefficients))] <- colnames(coefficients)[i]
  }
  # specify order of the levels to ensure that the coefficient output is in the same order as their input
  group <- factor(group, levels = colnames(coefficients))
  
  # Study indicator
  study <- rep(NA, length(yi))
  
  for (i in 1:dim(vcov)[3]) {
    study[(((i-1)*ncol(coefficients))+1):(i*ncol(coefficients))] <- rep(i,ncol(coefficients))
  }
  
  rma_dat <- data.frame(yi = yi, group = group, study = study)
  
  # Build the block matrix of within-study variances
  ws_var <- list()
  for (i in 1:dim(vcov)[3]) {
    ws_var[[i]] <- vcov[,,i]
  }
  S <- blockMatrixDiagonal(ws_var)
  
  ma.fit <- rma.mv(yi = yi, V = S, mods = ~ -1+group, random = ~ group|study, data = rma_dat,
                   struct= "UN", method = meta.method)
  
  out <- list(coefficients = ma.fit$beta[,1], # [,1] to coerce to vector with names.
              variances = diag(ma.fit$vb), # The estimated variances of the fixed effects
              se = sqrt(diag(ma.fit$vb)), # The estimated SEs of the fixed effects
              vcov = ma.fit$vb,
              psi = ma.fit$tau2)
  return(out)
}

which.abs.min <- function(x) 
  which.min(abs(x))

which.abs.max <- function(x)
  which.max(abs(x))

which.1 <- function(x)
  which.abs.min(x - 1)

#' @author Valentijn de Jong
#' @method unlist   listofperf
#' @export
unlist.listofperf <- function(x, ...) 
  sapply(x, `[[`, 1)

# Safe way of getting a function. 
# For some reason, this function can find functions that match.fun cannot.
# x function or character name thereof
# Returns function or error.
get.function <- function(x, ...) {
  if (is.function(x))
    return(x)
  else
    return(get(as.character(x), mode = "function"))
}

# Convert factor to binary
# In case the factor has more than 2 levels, by default the same occurs as in 
# glm: the first level is assumed to be the failure, and all others successes.
factor_as_binary <- function(x, failure_level = 1) {
  if (!is.factor(x)) stop(paste0("x must be a factor, x was ", class(x)))
  y <- rep(1, length(x))
  y[x == levels(x)[failure_level]] <- 0
  y[is.na(x)] <- NA
  y
}


