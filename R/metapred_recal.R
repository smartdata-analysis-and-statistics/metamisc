# computeRecal is an internal function for recalibrate() and metapred() to recalibrate the 
# intercept and/or coefficients.
# recalibrate() is to be used by predict.metapred() and is exported.
# object Model fit object, of class glm, lm or metapred
# newdata Data to use for recalibration.
# wholefit should whole fit be returned?

# To be replaced with one or two more basic functions: 1 for whole fit, 1 for returning coefficients only.

computeRecal <- function(object, newdata, b = NULL, f = ~ 1, estFUN = NULL, f.orig = NULL, wholefit = FALSE,  ...) {
  dots <- list(...)
  if (is.null(b)) b <- coef(object)
  if (is.null(estFUN)) {
    if (inherits(object, "metapred"))
      estFUN <- object$FUN$estFUN
    else estFUN <- as.character(class(object)[[1]])
  }
  estFUN <- match.fun(estFUN)
  
  # Convert outcome to binary, if it is a factor
  model_formula <- if (is.null(f.orig)) formula(object) else as.formula(f.orig)
  if (is.factor(newdata[, f2o(model_formula)]))
    newdata[ , f2o(model_formula)] <- factor_as_binary(newdata[ , f2o(model_formula)])
  
  # Make offset (linear predictor)
  X <- model.matrix(model_formula, data = newdata)
  lp <- X %*% b
  
  # offset must be in newdata.
  osdata   <- cbind(newdata, lp)
  f <- update.formula(formula(object), formula(f))
  if (is.null(object$family)) { 
    if (!is.null(dots[["family"]])) {
      refit <- estFUN(f, data = osdata, offset = lp, family = dots[["family"]])
    } else 
      refit <- estFUN(f, data = osdata, offset = lp)
  } else 
    refit <- estFUN(f, data = osdata, offset = lp, family = object$family)
  
  if (isTRUE(wholefit))
    return(refit)
  else
    return(coef(refit))
}

# p vector of predicted probs or outcomes
# y outcome vector
# estFUN estimation function
# family family
# which is "intercept" or "slope" or "add.slope".
# ... For compatibility only
pred.recal <- function(p, y, estFUN, family = gaussian, which = "intercept", ...) {
  if (is.character(family))  # Taken directly from glm()
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  estFUN <- match.fun(estFUN)
  lp <- family$linkfun(p)
  data <- data.frame(y = y, lp = lp)
  
  if (identical(which, "intercept"))
    out <- estFUN(formula = y ~ 1,  data = data, family = family, offset = lp)
  if (identical(which, "slope"))
    out <- estFUN(formula = y ~ lp + 1, data = data, family = family)
  if (identical(which, "add.slope"))
    out <- estFUN(formula = y ~ lp + 1, data = data, family = family, offset = lp)
  
  class(out) <- c("recal", class(out))
  out
}

# Object recal object
# parm "estimate" is the only viable obtion
# level 
# #' @export
# confint.recal <- function(object, parm = "estimate", level = .95, ...) { # Works for lm. And glm?
#   ses <- se(object, ...)
#   coefs <- coef(object, ...)
#   if(level < 0 || level > 1)
#     stop("Impossible confidence level. Possible levels: 0 < level < 1")
#   z <- qt(1 - (1 - level)/2, df = object$df.residual) # z = t distributed
#   data.frame("ci.lb" = coefs - z * ses, "ci.ub" = coefs + z * ses)
# }





## Split this into two functions.
# shrink <- function(object, newdata, method = "chisq", b = NULL, estFUN = NULL, ...) {
#   call <- match.call()
#   if (is.null(b)) b <- coef(object)
#   if (is.null(object$orig.coef))
#     object$orig.coef <- list()
#   if (!is.list(object$orig.coef))
#     stop("object is incompatible")
#   
#   if (identical(method, "chisq")) {
#     cs <- object$null.deviance - object$deviance
#     df <- object$df.null - object$df.residual
#     object$shrinkage.factor <- (cs - df) / cs
#   } 
#   else {
#   f <- formula( ~ lp)
#   
#   object$orig.coef[[length(object$orig.coef) + 1]] <- coef(object)
#   br <- computeRecal(object = object, newdata = newdata, f = f, estFUN = estFUN, ...)
#   
#   object$shrinkage.factor <- br[2] + 1
#   b[1] <- br[1]
#   b[-1] <- b[-1] * object$shrinkage.factor
#   object$coefficients <- b
#   }
# 
#   
#   if (is.call(object$call))
#   {
#     object$original.call <- object$call
#     object$call <- call
#   }
#   object
# }





# computeRecal(g, d3, estFUN = glm)

# \code{recalibrate} assumes coefficients are stored in \code{object$coefficients}
# and that \code{estFUN} accepts an \code{offset} argument.

#' Recalibrate a Prediction Model
#'
#' \code{recalibrate} is used to recalibrate a prediction model of classes \code{metapred, glm} or  \code{lm}.
#'
#' @param object A model fit object to be recalibrated, of class \code{metapred, glm} or \code{lm}, and more.
#' @param newdata data.frame containing new data set for updating.
#' @param f formula. Which coefficients of the model should be updated? Default: intercept only. Left-hand side may
#' be left out. See \link[stats]{formula} for details.
#' @param estFUN Function for model estimation. If left \code{NULL}, the function is automatically retrieved
#' for \code{metapred} objects. For other objects, the function with name corresponding to the first class
#' of the object is taken. E.g. \code{glm()} for \code{glm} objects.
#' @param ... Optional arguments to pass to \code{estFUN}.
#'
#' @details Currently only the coefficients are updated and the variances and other aspects are left untouched. 
#' For updating the entire model and all its statistics, see \link[stats]{update}.
#'
#' @return Recalibrated model fit object, of the same class as \code{object}. Generally, updated coefficients can
#' be retrieved with \code{coef()}.
#' 
#' @examples
#' data(DVTipd)
#' DVTipd$cluster <- 1:4 # Add a fictional clustering to the data set.
#' # Suppose we estimated the model in three studies: 
#' DVTipd123 <- DVTipd[DVTipd$cluster <= 3, ]
#' mp <- metamisc:::metapred(DVTipd123, strata = "cluster", f = dvt ~ vein + malign, 
#' family = binomial)
#' # and now want to recalibrate it for the fourth:
#' DVTipd4 <- DVTipd[DVTipd$cluster == 4, ]
#' metamisc:::recalibrate(mp, newdata = DVTipd4)
#'
#' @export
recalibrate <- function(object, newdata, f = ~ 1, estFUN = NULL, ...) {
  call <- match.call()
  if (is.null(object$orig.coef))
    object$orig.coef <- list()
  if (!is.list(object$orig.coef))
    stop("object is incompatible with recalibrate.")
  f <- as.formula(f)
  
  object$orig.coef[[length(object$orig.coef) + 1]] <- coef(object)
  br <- computeRecal(object = object, newdata = newdata, f = f, estFUN = estFUN, ...)
  i <- match(names(br), names(object$coefficients))
  object$coefficients[i] <- object$coefficients[i] + br
  
  if (is.call(object$call))
  {
    object$original.call <- object$call
    object$call <- call
  }
  object
}