#' Stacked Regression
#' 
#' This function combines one or more existing prediction models into a so/called meta-model.
#' 
#' @param models a list containing the historical prediction models, which can be defined in several ways. For instance, 
#' historical regression models can be specified using a named vector containing the regression coefficients of the 
#' individual predictors (no need to include the intercept term). List items may also represent an object for 
#' which the function \code{predict()} exists. 
#' @param family a description of the error distribution and link function to be used in the meta-model. This can be a character 
#' string naming a family function, a family function or the result of a call to a family function. (See \link[stats]{family} for 
#' details of family functions.)
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) 
#' containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, 
#' typically the environment from which \code{stackedglm} is called.
#' 
#' @keywords meta-analysis regression updating
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
#' 
stackedglm <- function(models, family = binomial, data) {
  call <- match.call()
  
  if (missing(models)) {
    stop("No historical models defined!")
  }
  if (!"list" %in% class(models)) {
    stop ("'models' should represent a list!")
  }
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  # Generate linear predictor of each model
  for (model in models) {
    if("numeric" %in% class(model)) {
      
    } else {
      # Models with a link function
      #try (predict(model, newdata=data, type="link"))
    }
  }

  
}

