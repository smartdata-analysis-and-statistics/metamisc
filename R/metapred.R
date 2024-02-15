# TODO 1: Thomas' ideas
# 2019-2-2: It would be nice to store some default performance estimates resulting from IECV in a "mm_perf" object (as as done in
#           'ccalc()' and 'oecalc'). It would then be easy to plot IECV results for these performance measures. 
#           Some suggested measures are: calibration slope, calibration-in-the-large and c-statistic. These performance
#           measures are fairly generic and should work for most glm-type models, including survival ones.
# @Valentijn: perhaps better to combine param perfFUN, genFUN and selFUN into one parameter with 3 or 4 distinct options
# @Valentijn: I suggest to omit intercept recalibration. The intercept issue can be addressed directly by specifying distinct 
# error functions that do or do not account for mis-calibration in intercept term

# TODO 2: metapred objects
# Add subset.metapred functionality for one.stage models.
# perfFUN = "cal.int" does not work for estFUN = "logistf", as logistf cannot take intercept only models
#     possible fix: use logistfirth intstead, and fix intercept-only models for that function.

# Bug when no family is specified for calibration slope:
# mp.slo <- metapred(DVTipd.reordered, strata = "cluster", perfFUN = "cal.slope", max.steps = 0)
# Error in pred.recal(p = p, y = y, estFUN = estFUN, family = family, which = "slope") : 
#   argument "family" is missing, with no default 

# add more output for metapred
# tol, aka tolerance for stopping.
# metapred.summaries docs (see glm.summaries {stats}), including subset
# add stepnumbers to print.metapred or print.fit
# Write tests for fitted.whatever
# add docs for fitted.whatever.
# check what happens for recalibrated mp.cv in fitted.mp.cv
# add rownames to data, to prevent error/mismatch in fitted.metapred/fitted.mp.cv # Maybe already fixed by as.data.frame()
# add the unchanged model to the next step's model list. Model comparison will be easier. Subsetting will be easier: this allows
#   a 'chain' of "changed" predictors to be added, thereby making the output clearer. 
# change "changed" for global models as well.
# add a new stratified.fit class, which contains mp.stratified.fit (or vice versa), such that subset(fit,
# select = "stratified") can work as intended. Right now, mp.stratified.fit  is a list that may only contain
# mp.stratum.fit objects.

# TODO 3: General
# add function: areTRUE. st.i == cl may otherwise lead to problems, if st.i has NA. Replaced with %in%. Now solved?
# variances for intercept recalibration.
# perf method
# penalties
# performance measurement, current is Brier for binomial
# add more options to penalty step
# Write test for new centerCovs function.

# TODO 4: Maybe later
# One-stage MA: predStep: is  type = response correct?
# centering: now only centers numeric variables. Should it also work for dummies of categorical variables?

# TODO 5: DONE:
# recal.int + interaction gives error. FIXED.
# add fitted.metapred(). DONE.
# categorical variables Done!
# Add summary.metapred. DONE.
# Added is.metapred()
# Added more robust tests for stepwise.
# Automatically remove observations with missing values.
# @Valentijn: I would center covariates by default, this often helps to improve generalizability and also speeds estimation
#   Valentijn: Done.
# @Valentijn: You can change default of meta.method to DL. This is a lot faster and has limited implications on estimated means.
#   Valentijn: Done.
# , cl.name in modelStep
# Remove Reduce() from perfStep, and make perfStep() compatible with multiple data sets (for cvFUN = fixed, or cvFUN = bs)
# predict.metapred: allow to generate predictions if y is not specified. 
#                   Currently, an error it thrown if setting dvt=NA in the example. This error traces back to the use 
#                   of model.matrix, which calls the function model.frame() that fails.
# Fixed by removing the outcome from the formula in the predict method.
# If max.steps and scope = formula are supplied, no stopping reason is given. SHould be fixed now (3 april 2019)

# TODO 6:

# predict.metapred: Allow so specify study for generating study-specific predictions 
# (as our IECV approach will yield a study-specific intercept term and regression coefficients)
# If no study is specified in the data frame, then predictions are based on the pooled coefficients.
# Change listwise deletion of missing data: add options to turn on/off, only remove observations with missing data within
# formula / scope.. 
# metapred currently cannot handle -1 in scope I think.

###### Outline
### Top: high-level functions
#
# metapred
# 
# mp.fit
# 
# mp.step + mp.global
# 
# mp.cv
# 
# mp.cv.val
# mp.cv.dev	+ mp.cv.recal
# 
# mp.cv.meta.fit + mp.recal.meta.fit
# 
# mp.meta.fit
#
# mp.stratified.fit
#
# mp.stratum.fit
#
### bottom : low level functions
###### 

#' Generalized Stepwise Regression for Prediction Models in Clustered Data
#'
#' Generalized stepwise regression for obtaining a prediction model that is validated 
#' with (stepwise) internal-external cross-validation, in or to obtain adequate 
#' performance across data sets. Requires data from individuals in multiple studies.
#' 
#' @author Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>
#' 
#' @references Debray TPA, Moons KGM, Ahmed I, Koffijberg H, Riley RD. A framework for developing, implementing, 
#' and evaluating clinical prediction models in an individual participant data meta-analysis. 
#' \emph{Stat Med}. 2013;32(18):3158-80. 
#'
#' @param data data.frame containing the data. Note that \code{metapred} removes observations with missing data \emph{listwise}
#' for all variables in \code{formula} and \code{scope}, to ensure that the same data is used in each model in each step. The
#' outcome variable should be \code{numeric} or coercible to numeric by as.numeric().
#' @param strata Character to specify the name of the strata (e.g. studies or clusters) variable
#' @param formula \code{formula} of the first model to be evaluated. \code{metapred} will start at \code{formula} and update it
#' using terms of \code{scope}. Defaults to full main effects model, where the first column in \code{data} is assumed to be
#' the outcome and all remaining columns (except \code{strata}) predictors. See \link[stats]{formula} for formulas in general.
#' @param estFUN Function for estimating the model in the first stage. Currently "lm", "glm" and "logistfirth" are supported.
#' @param scope \code{formula}. The difference between \code{formula} and \code{scope} defines the range of models examined in the 
#' stepwise search. Defaults to NULL, which leads to the intercept-only model. If \code{scope} is not nested in \code{formula}, 
#' this implies backwards selection will be applied (default). If \code{scope} is nested in \code{formula}, this implies forward 
#' selection will be applied. If equal, no stepwise selection is applied. 
#' @param retest Logical. Should added (removed) terms be retested for removal (addition)? \code{TRUE} implies bi-directional 
#' stepwise search.
#' @param max.steps Integer. Maximum number of steps (additions or removals of terms) to take. Defaults to 1000, which is
#' essentially as many as it takes. 0 implies no stepwise selection.
#' @param center logical. Should numeric predictors be centered around the cluster mean?
#' @param recal.int Logical. Should the intercept be recalibrated in each validation?
#' @param cvFUN Cross-validation method, on the study (i.e. cluster or stratum) level. "l1o" for leave-one-out cross-validation 
#' (default). "bootstrap" for bootstrap. Or "fixed", for one or more data sets which are only used for validation. A user written 
#' function may be supplied as well.
#' @param cv.k Parameter for cvFUN. For \code{cvFUN="bootstrap"}, this is the number of bootstraps. For \code{cvFUN="fixed"}, 
#' this is a vector of the indices of the (sorted) data sets. Not used for \code{cvFUN="l1o"}.
#' @param metaFUN Function for computing the meta-analytic coefficient estimates in two-stage MA. 
#' By default, \link[metafor]{rma.uni}, from the metafor package is used. Default settings are univariate random effects, 
#' estimated with "DL". Method can be passed trough the \code{meta.method} argument.
#' @param meta.method Name of method for meta-analysis. Default is "DL". For more options see \link[metafor]{rma.uni}.
#' @param predFUN Function for predicting new values. Defaults to the predicted probability of the outcome, using the link 
#' function of \code{glm()} or \code{lm()}.
#' @param perfFUN Function for computing the performance of the prediction models. Default: mean squared error 
#' (\code{perfFUN="mse"}).Other options are \code{"var.e"} (variance of prediction error), \code{"auc"} (area under the curve),
#' \code{"cal.int"} (calibration intercept), and \code{"cal.slope"} (multiplicative calibration slope) and \code{"cal.add.slope"}
#' (additive calibration slope).
#' @param genFUN Function or \code{list} of named functions for computing generalizability of the performance. 
#' Default: (absolute) mean (\code{genFUN="abs.mean"}). Choose \code{coef.var} for the coefficient of variation. If a \code{list},
#' only the first is used for model selection.
#' @param selFUN Function for selecting the best method. Default: lowest value for \code{genFUN}. Should be set to
#' "which.max" if high values for \code{genFUN} indicate a good model.
#' @param ... To pass arguments to estFUN (e.g. family = "binomial"), or to other FUNctions.
#'
#' @return A list of class \code{metapred}, containing the final model in \code{global.model}, and the stepwise
#' tree of estimates of the coefficients, performance measures, generalizability measures in \code{stepwise}.
#' 
#' @details Use \link{subset.metapred} to obtain an individual prediction model from a \code{metapred} object.
#' 
#'  Note that \code{formula.changes} is currently unordered; it does not represent the order of changes in the stepwise 
#'  procedure.
#'  
#'  \code{metapred} is still under development, use with care.
#' 
#' @examples 
#' data(DVTipd)
#' 
#'\dontrun{
#' # Explore heterogeneity in intercept and assocation of 'ddimdich'
#' glmer(dvt ~ 0 + cluster + (ddimdich|study), family = binomial(), data = DVTipd)
#'}
#' 
#' # Scope
#' f <- dvt ~ histdvt + ddimdich + sex + notraum
#' 
#' # Internal-external cross-validation of a pre-specified model 'f'
#' fit <- metapred(DVTipd, strata = "study", formula = f, scope = f, family = binomial)
#' fit
#' 
#' # Let's try to simplify model 'f' in order to improve its external validity
#' metapred(DVTipd, strata = "study", formula = f, family = binomial)
#' 
#' # We can also try to build a generalizable model from scratch
#' 
#'\dontrun{
#' # Some additional examples:
#' metapred(DVTipd, strata = "study", formula = dvt ~ 1, scope = f, family = binomial) # Forwards
#' metapred(DVTipd, strata = "study", formula = f, scope = f, family = binomial) # no selection
#' metapred(DVTipd, strata = "study", formula = f, max.steps = 0, family = binomial) # no selection
#' metapred(DVTipd, strata = "study", formula = f, recal.int = TRUE, family = binomial)
#' metapred(DVTipd, strata = "study", formula = f, meta.method = "REML", family = binomial)
#'}
#' # By default, metapred assumes the first column is the outcome.
#' newdat <- data.frame(dvt=0, histdvt=0, ddimdich=0, sex=1, notraum=0)
#' fitted <- predict(fit, newdata = newdat)
#' 
#' @references 
#' de Jong VMT, Moons KGM, Eijkemans MJC, Riley RD, Debray TPA. Developing more generalizable prediction models from pooled 
#' studies and large clustered data sets. \emph{Stat Med}. 2021;40(15):3533--59.
#' 
#' Riley RD, Tierney JF, Stewart LA. Individual participant data meta-analysis: a handbook for healthcare research. 
#' Hoboken, NJ: Wiley; 2021. ISBN: 978-1-119-33372-2.
#' 
#' Schmid CH, Stijnen T, White IR. Handbook of meta-analysis. First edition. Boca Raton: Taylor and Francis; 2020. ISBN: 978-1-315-11940-3.
#'   
#' @seealso  \code{\link{forest.metapred}}  for generating a forest plot of prediction model performance
#' @import stats
#'
#' @importFrom stats formula var
#'
#' @export
metapred <- function(data, strata, formula, estFUN = "glm", scope = NULL, retest = FALSE, max.steps = 1000, 
                     center = FALSE, recal.int = FALSE, cvFUN = NULL, cv.k = NULL,  # tol = 0,
                     metaFUN = NULL, meta.method = NULL, predFUN = NULL, perfFUN = NULL, genFUN = NULL,
                     selFUN = "which.min",
                     ...) {
  call <- match.call()
  
  # Formula
  if (missing(formula) || is.null(formula)) formula <- stats::formula(data[ , -which(colnames(data) == strata)])  
  if (is.null(scope)) scope <- f2iof(formula)
  updates <- getFormulaDiffAsChar(formula, scope)
  
  # Data
  data <- get_all_vars(formula = addg2f(formula, scope, terms = strata), data = data) # drop unnecessary vars
  data <- droplevels(remove.na.obs(as.data.frame(data)))  # drop observations with missings in remaining vars
  
  if (is.factor(data[, f2o(formula)]))
    data[ , f2o(formula)] <- factor_as_binary(data[ , f2o(formula)])
  
  # One- vs two-stage, stratification and centering
  if(is.null(two.stage <- list(...)$two.stage) ) two.stage <- TRUE
  strata.i <- as.vector(data[, strata])
  strata.u <- sort(unique(strata.i))
  if (center)
    data <- centerCovs(data = data, y.name = f2o(formula), cluster.name = strata)
  
  # Functions
  if (is.null(cvFUN))   cvFUN   <- l1o
  if (is.null(metaFUN)) metaFUN <- urma
  if (is.null(perfFUN)) perfFUN <- "mse"
  if (is.null(genFUN))  genFUN  <- abs.mean
  if (is.null(meta.method)) meta.method <- "DL"
  # Change to "-" when perfFUN <- R2 or some other measure for which greater = better.
  
  estFUN.name <- estFUN
  estFUN  <- get.function(estFUN)
  cvFUN   <- get.function(cvFUN)
  # perfFUN <- get(perfFUN)
  # genFUN  <- get(genFUN) # now happens in mp.cv.val
  selFUN  <- get.function(selFUN)
  metaFUN <- get.function(metaFUN)
  
  # genFUN.add <- dots[["genFUN.add"]] 
  # dots[["genFUN.add"]] <- NULL
  
  # Folds
  folds <- cvFUN(strata.u, k = cv.k)
  if (!isTRUE(length(folds[["dev"]]) > 0) || !isTRUE(length(folds[["dev"]][[1]]) > 0))
    stop("At least 1 cluster must be used for development.")
  
  # Fitting
  fit <- mp.fit(formula = formula, data = data, remaining.changes = updates, st.i = strata.i, st.u = strata.u, folds = folds,
                recal.int = recal.int, retest = retest, max.steps = max.steps, tol = 0,
                estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, predFUN = predFUN, perfFUN = perfFUN,
                genFUN = genFUN, selFUN = selFUN, ...)
  
  # mp.args <- c(list(formula = formula, data = data, remaining.changes = updates, st.i = strata.i, st.u = strata.u, folds = folds,
  #                   recal.int = recal.int, retest = retest, max.steps = max.steps, tol = 0, 
  #                   estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, predFUN = predFUN, perfFUN = perfFUN,
  #                   genFUN = genFUN, selFUN = selFUN), dots)
  # 
  # fit <- do.call(mp.fit, args = mp.args ) 
  
  predFUN <- getPredictMethod(fit$stepwise$s0$cv$unchanged, two.stage = two.stage, predFUN = predFUN)
  formula.final <- fit$global.model$formula
  
  out <- c(fit, list(call = call, strata = strata, data = data, folds = folds, # add nobs and strata.nobs
                     formula.start = formula, scope = scope, formula = formula.final,
                     formula.changes = getFormulaDiffAsChar(formula.final, formula), 
                     # NOTE: formula.changes is currently unordered!
                     options = list(cv.k = cv.k, meta.method = meta.method, recal.int = recal.int,
                                    center = center, max.steps = max.steps, retest = retest, two.stage = two.stage), # add: tol
                     FUN = list(cvFUN = cvFUN, predFUN = predFUN, perfFUN = get.function(perfFUN), metaFUN = metaFUN, genFUN = genFUN, 
                                selFUN = selFUN, estFUN = estFUN, estFUN.name = estFUN.name)))
  class(out) <- c("metapred")
  return(out)
}

# #' The function \code{summary} can be used to obtain an extensive summary of the stepwise process and the final model. 
# #' \code{subset(x, select = "best.cv")} can be used to obtain a summary of the model with best generalizability in the
# #' cross-validation, whereas \code{subset(x, select = "global")} gives a summary of the global model (i.e. fitted on all 
# #' strata).

# For prediction of newdata. 
# Object metapred object
# newdata data.frame, defaults to NULL, which implies the fitted data
# strata character, name of strata variable in newdata. Defaults to name in fitted object.
# type character. Type of prediction. This is intended to override the default of glm and lm.
# ("response" or "link" possible; "terms" not implemented)
# recal.int logical. Recalibrate the intercept before prediction? Defaults to same as development for metapred object,
# center logical, Center covariates, before prediction? Defaults to same as development of metapred object,
# ... For compatibility only.
#' @author Valentijn de Jong
#' @importFrom stats predict
#' @method predict   metapred
#' @export
predict.metapred <- function(object, newdata = NULL, strata = NULL, type = "response", 
                             recal.int = NULL, center = NULL, ...) {
  if (is.null(newdata))
    newdata <- object[["data"]]
  if (is.null(strata))
    strata <- object$strata
  if (is.null(center))
    center <- object$options$center
  if (center)
    newdata <- centerCovs(data = newdata, y.name = f2o(formula(object)), cluster.name = strata)
  if (is.null(recal.int))
    recal.int <- object$options$recal.int
  if (isTRUE(recal.int))
    object <- recalibrate(object = object, newdata = newdata)
  
  object$FUN$predFUN(object = object, newdata = newdata, type = type, ...)[, 1] # [, 1] such that a vector is returned.
}

#' Extract Model Fitted Values
#' 
#' Extract the fitted values of a \code{metapred} object. By default returns fitted values of the model in the 
#' cross-validation procedure.
#' 
#' Function still under development, use with caution.
#' 
#' Only returns type = "response".
#' 
#' @author Valentijn de Jong
#' @importFrom stats fitted
#' @method fitted   metapred
#' @param object object of class metapred 
#' @param select character. Select fitted values from "cv" (default) or from "global" model.
#' @param step character or numeric. Name or number of step to select if \code{select} = "cv". Defaults to best step.
#' @param model character or numeric. Name or number of model to select if \code{select} = "cv". Defaults to
#' best model.
#' @param as.stratified logical. \code{select} = "cv" determines whether returned predictions are stratified in a list 
#' (\code{TRUE}, default) or in their original order (\code{FALSE}).
#' @param type character. Type of fitted value.
#' @param ... For compatibility only.
#' @export
fitted.metapred <- function(object, select = "cv", step = NULL, model = NULL, 
                            as.stratified = TRUE, type = "response", ...) {
  if (isTRUE(select == "cv")) {
    ftd <- fitted(subset.metapred(x = object, select = select, step = step, model = model, type = type, ...))
    if (as.stratified)
      return(ftd)
    ftd.v <- Reduce(rbind, ftd) #as vector
    return(ftd.v[match(rownames(ftd.v), rownames(object[["data"]])) ] ) # return to original ordering.
  }
  
  if (isTRUE(select == "global"))
    return(predict.metapred(object = object, newdata = NULL, type = type, ...))
  stop("select must equal 'cv' or 'global'.")
}

# #' Only returns type = "response".
# #' @author Valentijn de Jong
# #' @importFrom stats residuals
# #' @method residuals   metapred
# #' @export
# residuals.metapred <- function(object, select = "cv", step = NULL, model = NULL, as.stratified = TRUE, ...) {
#   y <- object$data[ , f2o(formula(object))]
#   ftd <- fitted.metapred(object = object, select = select, step = step, 
#                        model = model, as.stratified = as.stratified, ...)
#   if (as.stratified) {
#     ftd <- fitted(mp, as.stratified = TRUE)
#     y <- mp$data[ , f2o(formula(mp))]
#     
#     st.i <- mp[["data"]][[mp[["strata"]]]]
#     # now somehow sort them into a list with similar dim and dimnames as ftd, using st.i
#     
#     stop("To be implemented")
#   } else {
#     y - ftd
#   }
# }

# #' Extract the regression coefficients
# #' 
# #' The \code{coef} function extracts the estimated coefficients of the final model of objects of class \code{"metapred"}.
# #' @param object A fitted \code{metapred} object
# #' @param \ldots Optional arguments (currently not supported).
# #' 
# #' @method coef metapred
# #' @export
#' @author Valentijn de Jong
#' @method coef   metapred
#' @importFrom stats coef
#' @export
coef.metapred <- function(object, stratified = FALSE, ...) {
  if (stratified) 
    return(coef(object$global.model$stratified.fit))
  return(coef(object$global.model))
}

# Note returned formula is for the final model!
#' @author Valentijn de Jong
#' @importFrom stats formula
#' @method formula   metapred
#' @export
formula.metapred <- function(x, ...)
  x$global.model$formula

# Returns family or NULL
#' @author Valentijn de Jong
#' @importFrom stats family
#' @method family   metapred
#' @export
family.metapred <- function(object, ...) {
  if (!is.null((f <- object$stepwise$s0$cv[[1]]$family)))
    f
  else if (!is.null(f <- object$family))
    f
  else
    NULL
}

#' @author Valentijn de Jong
#' @method print   metapred
#' @export
print.metapred <- function(x, ...) {
  cat("Call: ")                       
  print(x$call) # cat cannot handle a call
  cat("\n") 
  print.mp.fit(x)
}

#' @author Valentijn de Jong
#' @method summary   metapred
#' @export
summary.metapred <- function(object, ...) {
  cat("Call: ")                       
  print(object$call) # cat cannot handle a call
  cat("\n") 
  summary.mp.fit(object)
}

# Test whether object is meta.pred
#' @author Valentijn de Jong
#' @importFrom methods is
#' @method is   metapred
#' @export
is.metapred <- function(object)
  inherits(object, "metapred")

# Implementation of the subset method
# #' @author Valentijn de Jong
# #' @method subset   metapred
# #' @export
# old.subset.metapred <- function(x, select = "best.cv", ...) {
#   if (identical(select, "best.cv")) 
#     return(mp.step.get.best(x[["stepwise"]][[x[["best.step"]]]]))
#   if (identical(select, "global"))
#     return(x[["global.model"]])
#   stop("select must equal 'best.cv' or 'global'.")
# }

#' Subsetting metapred fits
#' 
#' Return a model from the cross-validation procedure or the final 'global' model. Caution: This function is 
#' still under development.
#' 
#' @author Valentijn de Jong
#' 
#' @param x metapred object
#' @param select Which type of model to select: "cv" (default), "global", or (experimental) "stratified", 
#' or "stratum".
#' @param step  Which step should be selected? Defaults to the best step. 
#' numeric is converted to name of the step: 0 for an unchanged model, 1 for the first change... 
#' @param model Which model change should be selected? NULL (default, best change) or character name of variable
#' or (integer) index of model change. 
#' @param stratum Experimental. Stratum to return if select = "stratum".
#' @param add Logical. Add data, options and functions to the resulting object? Defaults to \code{TRUE}. 
#' Experimental.
#' @param ... For compatibility only.
#' 
#' @examples 
#' data(DVTipd)
#' DVTipd$cluster <- letters[1:4] # Add a fictional clustering to the data.
#' mp <- metapred(DVTipd, strata = "cluster", formula = dvt ~ histdvt + ddimdich, family = binomial)
#' subset(mp) # best cross-validated model
#' subset(mp, select = "global") # Final model fitted on all strata.
#' subset(mp, step = 1) # The best model of step 1
#' subset(mp, step = 1, model = "histdvt") # The model in which histdvt was removed, in step 1.
#' 
#' @return An object of class \code{mp.cv} for select = "cv" and an object of class \code{mp.global} for select = "global". 
#' In both cases, additional data is added to the resulting object, thereby making it suitable for further methods.
#' @export
subset.metapred <- function(x, select = "cv", step = NULL, model = NULL, stratum = NULL, add = TRUE, ...) {
  # Global model
  if (identical(select, "global"))
    out <- x[["global.model"]]
  else {
    # Same for cv, stratified and stratum:
    # Step
    if (is.null(step))
      step <- x$best.step
    if (is.numeric(step))
      step <- getStepName(step)
    
    # Model
    if (is.null(model))
      model <- x[["stepwise"]][[step]][["best.change"]]
    
    # Select fit:
    if (identical(select, "cv") || identical(select, "iecv"))                        # cv / iecv fit
      out <- x[["stepwise"]][[step]][["cv"]][[model]]
    else if (identical(select, "stratified") || identical(select, "stratified.fit")) # stratified fit
      return(x[["stepwise"]][[step]][["cv"]][[model]][["stratified.fit"]])
    else if (identical(select, "stratum")    || identical(select, "strata"))         # stratum fit
      out <- x[["stepwise"]][[step]][["cv"]][[model]][["stratified.fit"]][[stratum]]
    else
      stop("select must equal 'cv', 'global', 'stratified' or 'stratum'.")
  }   
  
  # This is the real reason for why this function exists. To add this stuff to a mp.cv or mp.global object.
  # Normally those do not have this data, for memory/performance reasons.
  # In hindsight: This is actually a very useful function, as it makes finding a model much easier.
  if (add) {
    out$data    <- x$data
    out$strata  <- x$strata
    out$folds   <- x$folds
    out$FUN     <- x$FUN
    out$options <- x$options
  }
  
  out
}

# Perform fit a model using cross-validation, for metapred
# formula formula to start with
# data data.frame, containing dev and val data
# remaining.changes predictor terms to add or remove
# st.i numeric vector, cluster indicators
# st.u unique values in st.i
# folds list, folds as in utils
# recal.int logical, recalibrate intercept?

# retest logical, should removed (added) predictors be added (removed) again?
# max.steps numeric, maximum number of steps (predictor additions/removals) to be taken
# tol numeric, tolerance, minimum change in generalizability to accept a different model # To be implemented

# estFUN function, used for obtaining predict method
# predFUN function, user supplied predict method, overrides estFUN's predict()
# perfFUN function, function for computing performance, defaults to mse = mean square error
# genFUN function, function for computing generalizability, defaults to absolute mean.
# ... options for predictMethod, perfFUN, and genFUN.
# Returns object of class mp.cv.val, which is a validated mp.cv.dev
mp.fit <- function(formula, data, remaining.changes, st.i, st.u, folds, recal.int = FALSE, 
                   retest = FALSE, max.steps = 3, tol = 0,
                   estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                   perfFUN = mse, genFUN = abs.mean, selFUN = which.min,
                   two.stage = TRUE, ...) {
  out <- steps <- list()
  
  ## Step 0
  # As remaining.changes = c("") yields unchanged formula.
  step.count <- 0
  steps[[getStepName(step.count)]] <- mp.step(formula = formula, data = data, remaining.changes = c(""), 
                                              st.i = st.i, st.u = st.u, folds = folds, recal.int = recal.int, 
                                              retest = FALSE, two.stage = two.stage,
                                              estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, predFUN = predFUN, 
                                              perfFUN = perfFUN, genFUN = genFUN, selFUN = selFUN, ...)
  steps[[getStepName(step.count)]][["step.count"]] <- step.count
  out[["best.step"]] <- getStepName(step.count)
  out[["stop.reason"]] <- "no changes were requested."
  current.model <- mp.step.get.best(steps[[1]])
  current.model[["remaining.changes"]] <- remaining.changes
  gen.diff <- Inf # TBI
  
  if (!identical(length(remaining.changes), 0L))
    repeat {
      ## Loop management
      if (isTRUE(gen.diff <= tol)) { # TBI
        out[["stop.reason"]] <- "improvement <= tolerance."
        break
      }
      if (isTRUE(length(current.model[["remaining.changes"]]) <= 0)) {
        out[["stop.reason"]] <- "all changes were implemented."
        break
      }
      if (isTRUE(step.count >= max.steps)) {
        out[["stop.reason"]] <- "max.steps was reached."
        break
      }
      step.count <- step.count + 1
      
      ## Perform a step
      steps[[getStepName(step.count)]] <- mp.step(formula = current.model[["formula"]], data = data,
                                                  remaining.changes = current.model[["remaining.changes"]],
                                                  st.i = st.i, st.u = st.u, folds = folds, recal.int = recal.int,
                                                  retest = retest, two.stage = two.stage,
                                                  estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, predFUN = predFUN,
                                                  perfFUN = perfFUN, genFUN = genFUN, selFUN = selFUN, ...)
      steps[[getStepName(step.count)]][["step.count"]] <- step.count
      ## Model selection
      # This step
      best.new.model <- mp.step.get.best(steps[[getStepName(step.count)]])
      
      # Overall
      if (mp.which.best.change(list(current.model, best.new.model), selFUN = selFUN) <= 1) { # TBI: 0 for gen.diff < tol
        out[["stop.reason"]] <- "no improvement was possible."
        out[["best.step"]] <- getStepName(step.count - 1)
        break
      } else {
        current.model <- best.new.model
        out[["best.step"]] <- getStepName(step.count)
      }
    }
  
  # Return a global model and the stepwise list
  if (is.null(two.stage) || isTRUE(two.stage)) {
    out[["global.model"]] <- mp.global.2st(current.model, metaFUN = metaFUN, meta.method = meta.method)
    out[["nobs.strata"]] <- sapply(out[["global.model"]][["stratified.fit"]], nobs)
  } else {
    out[["global.model"]] <- mp.global.1st(formula = formula(current.model), data = data,
                                           estFUN = estFUN, ...) # strata = strata, to be added for random effects.
  }
  
  out[["stepwise"]] <- steps
  out[["step.count"]] <- step.count
  out[["nobs"]] <- sum(out[["nobs.strata"]])
  class(out) <- "mp.fit"
  
  out
}

#' @author Valentijn de Jong
#' @method print   mp.fit
#' @export
print.mp.fit <- function(x, ...) {
  for (i in seq_along(x[["stepwise"]])) {
    if (i == 1) {
      cat("Started with model:\n")
      print(x$stepwise[[i]][["start.formula"]])
    } else if (i > 2) {
      cat("\n")
      cat("Continued with model:\n")
      print(x$stepwise[[i]][["start.formula"]])
    }
    print.mp.step(x$stepwise[[i]], show.f = FALSE)
  }
  
  cat("\n")
  cat("Cross-validation stopped after", x[["step.count"]], "steps, as", x[["stop.reason"]])
  cat(" Final model:\n")
  print(x[["global.model"]])
}

#' @author Valentijn de Jong
#' @method summary   mp.fit
#' @export
summary.mp.fit <- function(object, ...) {
  for (s in seq_along(object[["stepwise"]]))
    summary.mp.step(object[["stepwise"]][[s]])
  cat("\n")
  cat("Cross-validation stopped after", object[["step.count"]], "steps, as", object[["stop.reason"]])
  cat(" Final model:\n")
  print(object[["global.model"]])
}

# "Method" for selecting best mp.cv from a mp.cv object
# step mp.step
# Returns mp.cv object, with best value for genFUN.
# mp.step.best.change <- function(step, ...)
#   step[["cv"]][[which(sapply(step[["cv"]], `[[`, "best.change"))]]

# Select best model
mp.which.best.change <- function(cvs, selFUN = which.min, ...)
  selFUN(sapply(cvs, `[[`, "gen"))

mp.step.get.best <- function(step, selFUN = which.min, ...) 
  step[["cv"]][[mp.which.best.change(step[["cv"]], selFUN = selFUN, ...)]]

mp.step.get.change <- function(step, ...)
  mp.step.get.best(step)[["changed"]]

# Perform one step in the fitting process of mp.fit
# formula formula to start with
# data data.frame, containing dev and val data
# remaining.changes predictor terms to add or remove
# st.i numeric vector, cluster indicators
# st.u unique values in st.i
# folds list, folds as in utils
# recal.int logical, recalibrate intercept?
# two.stage logical, is it a two stage model? For future use.
# estFUN function, used for obtaining predict method
# predFUN function, user supplied predict method, overrides estFUN
# perfFUN function, function for computing performance, defaults to mse = mean square error
# genFUN function, function for computing generalizability, defaults to absolute mean.
# ... options for predictMethod, perfFUN, and genFUN.
# Returns object of class mp.cv.val, which is a validated mp.cv.dev
mp.step <- function(formula, data, remaining.changes, st.i, st.u, folds, recal.int = FALSE, 
                    two.stage = TRUE, retest = FALSE,
                    estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                    perfFUN = mse, genFUN = abs.mean, selFUN = which.min, ...) {
  cv <- out <- list()
  out[["start.formula"]] <- formula
  
  for (fc in seq_along(remaining.changes) )
  {
    change <- remaining.changes[fc]
    
    # Produce formula for changes and no changes:
    if (identical(remaining.changes, "")) {
      name <- "unchanged"
      new.formula <- formula
    } else {
      name <- as.character(change)
      new.formula <- updateFormula(formula, change) 
    }
    
    # Run
    # change.name is now saved here and used elsewhere. To be simplified..
    cv[[name]] <- mp.cv(formula = new.formula, data = data, st.i = st.i, st.u = st.u,
                        folds = folds, recal.int = recal.int, two.stage = two.stage,
                        estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method,
                        predFUN = predFUN, perfFUN = perfFUN, genFUN = genFUN, change = change, ...) 
    # Save changes
    cv[[name]][["remaining.changes"]] <- if (retest) remaining.changes else remaining.changes[-fc]
    # cv[[name]][["changed"]] <- change
  }
  
  out[["best.change"]] <- mp.which.best.change(cv, selFUN = selFUN)
  out[["cv"]] <- cv
  
  class(out) <- "mp.step"
  out
}

#' @author Valentijn de Jong
#' @method summary   mp.step
#' @export
summary.mp.step <- function(object, ...) {
  cv <- object[["cv"]]
  
  # Model comparison is handled by print.mp.step()
  if (!is.null(object[["step.count"]]))
    cat("\n", "Step ", object[["step.count"]], ". ", sep = "")
  if (identical(cv[[1]][["changed"]], "")) {
    # cat("\n")
    cat("Tested 1 model in this step. ") 
  } else {
    # cat("\n")
    cat("Tested", length(cv), "models in this step. ")
    print.mp.step(object)
    cat("\n")
  }
  
  # Model coefficients are handled here, using print.mp.cv()
  for (i in seq_along(cv)) {
    ch <- cv[[i]][["changed"]]
    if (identical(ch, "")) ch <- "nothing"
    cat("Changing", ch, "yields:\n")
    print.mp.cv(cv[[i]])
    if (isTRUE(i == object[["best.change"]][[1]]))
      cat("This model has best generalizability in this step.")
    cat("\n")
  }
}

#' @author Valentijn de Jong
#' @method print   mp.step
#' @export
print.mp.step <- function(x, show.f = TRUE, ...) {
  if (show.f) {
    cat("Starting with model:\n")
    print(x[["start.formula"]])
  }
  cat("\n")
  cat("Generalizability:\n")
  print(data.frame(lapply(x[["cv"]], `[[`, "gen.all")))
}

# Turn a cross-validated model into a full or 'global' model
# mp.cv.dev object of class mp.cv.dev
# metaFUN function for meta-analysis
# meta.method option for metaFUN
# ... optional arguments for metaFUN
mp.global.2st <- function(cv.dev, metaFUN = urma, meta.method = "DL") {
  
  out <- c(cv.dev, mp.meta.fit(stratified.fit = cv.dev[["stratified.fit"]], 
                               metaFUN = metaFUN, meta.method = meta.method) )
  
  class(out) <- c("mp.global", class(out))
  out
}

mp.global.1st <- function(formula, data, strata = NULL, estFUN = "glm", ...) {
  # TODO for glmer / lmer: add random intercepts to formula here.
  dots <- list(...)
  dots$two.stage <- NULL
  args <- c(list(formula = formula, data = data), dots)
  fit <- do.call(estFUN, args)
  out <- mp.stratum.fit(fit)
  out$formula <- formula
  out[["nobs.strata"]] <- nobs(out) # To be removed??
  class(out) <- c("mp.global", "mp.stratum.fit", class(out)) # "mp.global.1st" ?! or not?
  out
}

# #' @author Valentijn de Jong
# #' @method family   mp.global
# #' @export
# family.mp.global <- function(object, ...)
#   object$family

# Fitted? Predict?

#' @author Valentijn de Jong
#' @method predict   mp.global
#' @export
predict.mp.global <- function(object, newdata = NULL, strata = NULL, type = "response", 
                              recal.int = NULL, center = NULL, ...) 
  predict.metapred(object = object, newdata = newdata, strata = strata, type = type, 
                   recal.int = recal.int, center = center, ...)

#' @author Valentijn de Jong
#' @method print   mp.global
#' @export
print.mp.global <- function(x, ...) {
  cat("Meta-analytic model of prediction models estimated in", x$n.clusters, "strata. Coefficients:" , "\n")
  print(x[["coefficients"]])
}

#' @author Valentijn de Jong
#' @method summary   mp.global
#' @export
summary.mp.global <- function(object, ...) {
  print.mp.global(object, ...)
  cat("\n")
  print.mp.cv(object, ...)
}

# Make a new meta-model and cross-validate it, for metapred
# formula formula
# data data.frame
# st.i vector, cluster indicators. char or numeric
# st.u unique values in st.i
# recal.int logical, recalibrate intercept?
# estFUN function, for estimating, e.g. glm
# metaFUN function, for producing meta model
# meta.method character, option for metaFUN
# predFUN function, user supplied predict method, overrides estFUN
# perfFUN function, function for computing performance, defaults to mse = mean square error
# genFUN function, function for computing generalizability, defaults to absolute mean.
# ... options for predictMethod, perfFUN, and genFUN.
# Returns object of classes mp.cv, mp.cv.val, mp.cv.dev, which is a list of meta-analytic models developed on dev folds,
# and a validated on val folds 
mp.cv <- function(formula, data, st.i, st.u, folds, recal.int = FALSE, two.stage = TRUE,
                  estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                  perfFUN = mse, genFUN = abs.mean, change = NULL, ...) {
  out <- mp.cv.dev(formula = formula, data = data, st.i = st.i, st.u = st.u, folds = folds, two.stage = two.stage,
                   estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, change = change, ...)
  
  out <- mp.cv.val(cv.dev = out, data = data, st.i = st.i, folds = folds, recal.int = recal.int, two.stage = two.stage,
                   estFUN = estFUN, predFUN = predFUN, perfFUN = perfFUN, genFUN = genFUN, ...)
  
  class(out) <- c("mp.cv", class(out))
  out
}

# Obtain fitted values of a (\code{subset}ted from metapred) mp.cv object
# object mp.cv
# returns list of fitted values, with length equal to number of strata. 
#' @author Valentijn de Jong
#' @method fitted   mp.cv
#' @export
fitted.mp.cv <- function(object, ...) {
  data <- object[["data"]]
  if (is.null(data))
    stop("Use subset(metapred()) to obtain the cv.")
  strata <- object$strata
  center <- object$options$center
  if (center)
    data <- centerCovs(data = data, y.name = f2o(formula(object)), cluster.name = strata)
  
  fitted.mp.cv.dev(object, data = data, folds = object$folds, st.i = data[[object$strata]], ...)
}

# fitted.mp.cv.val <- function(object, newdata, folds, st.i, recal.int = FALSE, predFUN = NULL, ...) 
#   fitted.mp.cv.dev(object = object, newdata = newdata, folds = folds, st.i = st.i, predFUN = predFUN)

#' @author Valentijn de Jong
#' @method print   mp.cv
#' @export
print.mp.cv <- function(x, ...) {
  print.mp.cv.dev(x, ...)
  cat("\n")
  print.mp.cv.val(x, ...)
}

# Validate a mp.cv.dev
# cv.dev, mp.cv.dev
# data data.frame, containing old and new data
# st.i numeric vector, cluster indicators
# folds list, folds as in utils
# recal.int logical, recalibrate intercept?
# two.stage logical, is it a two stage model?
# estFUN function, used for obtaining predict method
# predFUN function, user supplied predict method, overrides estFUN
# perfFUN function, function for calculating performance, defaults to mse = mean square error
# genFUN function, function for calculating generalizability, defaults to absolute mean.
# add.perfFUN list of functions, additional performance functions. Evaluated but not used for selection.
# add.genFUn list of functions, additional generalizability functions. Evaluated but not used for selection.
# ... options for predictMethod, perfFUN, and genFUN.
# Returns object of class mp.cv.val, which is a validated mp.cv.dev
mp.cv.val <- function(cv.dev, data, st.i, folds, recal.int = FALSE, two.stage = TRUE, 
                      estFUN = glm, predFUN = NULL, perfFUN = mse, 
                      genFUN = abs.mean, plot = F, ...) {
  dots <- list(...)
  pfn <- if (is.character(perfFUN)) perfFUN else "Performance"
  cv.dev[["perf.name"]] <- pfn # To be removed!??!!?
  # perfFUN <- match.fun(perfFUN)
  
  # Recalibrate?
  if (isTRUE(recal.int))
    cv.dev <- mp.cv.recal(cv.dev = cv.dev, newdata = data, estFUN = estFUN, folds = folds)
  
  # Predict outcome
  p <- fitted.mp.cv.dev(object = cv.dev, data = data, folds = folds, st.i = st.i, predFUN = predFUN, two.stage = two.stage)
  cv.dev[["nobs.val"]] <- sapply(p, length)
  
  # Necessary for performance computation
  outcome <- f2o(formula(cv.dev))
  
  # Performance
  if (!is.list(perfFUN)) 
    perfFUN <- list(perfFUN)
  
  # Names of Performance measures
  if (identical(length(names(perfFUN)), length(perfFUN))) {
    perf.names <- names(perfFUN)
  } else if (all(sapply(perfFUN, is.character))) {
    perf.names <- unlist(perfFUN)
  } else
    perf.names <- as.character(seq_along(perfFUN))
  cv.dev[["perf.names"]] <- perf.names
  
  # Multiple performance measures may be calculated.
  perfcalc <- function(perfFUN, cv.dev, folds, outcome, st.i, data, estFUN, p) {
    perfFUN <- match.fun(perfFUN)
    
    perf.full <- perf.str <- list()
    
    for (i in seq_len(cv.dev[["n.cv"]])) {
      perf.str[[length(perf.str) + 1]]   <- getclName(folds[["val"]][[i]])
      perf.full[[length(perf.full) + 1]] <- perfFUN(p[[i]], data[folds[["val"]][[i]] == st.i, outcome], data = data,
                                                    fit = cv.dev[["cv"]][[i]], estFUN = estFUN, ...)
    }
    names(perf.full) <- apply(as.matrix(folds[["val"]]), 1, getclName)
    
    out <- data.frame(val.strata = unlist(perf.str), 
                      estimate = unlist(sapply(perf.full, `[[`, 1)), 
                      se = NA, var = NA, ci.lb = NA, ci.ub = NA, 
                      measure = pfn,
                      n = unlist(cv.dev[["nobs.val"]]),
                      class = unlist(lapply(lapply(perf.full, class), '[[', 1)))
    tryCatch(
      out[["var"]] <- unlist(lapply(perf.full, variances)),  
      error = function(e) print(paste("Skipping variance estimation for", class(perf.full[[1]])[[1]], sep = " ")))
    tryCatch(
      out[["se"]] <- sqrt(out[["var"]]),
      error = function(e) print(paste("Skipping se estimation for", class(perf.full[[1]])[[1]], sep = " ")))
    tryCatch(
      out[, c("ci.lb", "ci.ub")] <- unlist(t(sapply(perf.full, get_confint))),  
      error = function(e) print(paste("Skipping ci estimation for", class(perf.full[[1]])[[1]], sep = " ")))
    
    row.names(out) <- names(cv.dev[["cv"]])
    out[, "estimate"]  <- as.numeric(out[, "estimate"]) 
    class(out) <- c("perf", class(out))
    out # cv.dev[["perf"]] = out
  }
  
  perf.all <- list()
  for (fun.id in seq_along(perfFUN)) # Single brackets intended!
    perf.all[[fun.id]] <- perfcalc(perfFUN[[fun.id]], cv.dev = cv.dev, folds = folds, outcome = outcome, 
                                   st.i = st.i, data = data, estFUN = estFUN, p = p)
  
  names(perf.all) <- perf.names
  cv.dev[["perf.all"]] <- perf.all   # Future compatibility  
  cv.dev[["perf"]] <- perf.all[[1]]  # Backwards compatibility
  
  # Generalizibility
  if (!is.list(genFUN)) 
    genFUN <- list(genFUN)
  
  # Names of generalizability measures
  if (identical(length(names(genFUN)), length(genFUN))) {
    gen.names <- names(genFUN)
  } else if (all(sapply(genFUN, is.character))) {
    gen.names <- unlist(genFUN)
  } else
    gen.names <- as.character(seq_along(genFUN))
  
  # Computation of generalizability
  gen.all <- rep(NA, length(genFUN))
  
  for (fun.id in seq_along(genFUN)) { # Single brackets intended!
    genfun <- match.fun(genFUN[[fun.id]])
    args <- c(list(object = cv.dev[["perf"]],
                   coef = coef(cv.dev[["stratified.fit"]]),
                   title = paste("Model change: ~", cv.dev[["changed"]]),
                   xlab = as.character(pfn)
    ), dots)
    if (!is.null(fam <- cv.dev$family))
      args$family <- fam
    if (two.stage)
      args$coef.se <- se(cv.dev[["stratified.fit"]]) # is compatible with two.stage = T, but not with = F
    
    gv <- do.call(genfun, args = args )
    
    # gv <- genfun(cv.dev[["perf"]], coef = coef(cv.dev[["stratified.fit"]]), coef.se = se(cv.dev[["stratified.fit"]]),
    # title = paste("Model change: ~", cv.dev[["changed"]]), xlab = as.character(pfn), ...)
    # As 'foo[[bar]] <- NULL' is not allowed # Also, the rest of the code does not expect a plot yet.
    gen.all[[fun.id]] <- if (is.null(gv) || !is.numeric(gv) || !identical(length(gv), 1L)) NaN else gv 
  }
  
  names(gen.all) <- gen.names
  cv.dev[["gen.all"]] <- gen.all
  cv.dev[["gen"]] <- gen.all[[1]]
  
  class(cv.dev) <- c("mp.cv.val", class(cv.dev))
  if (plot)
    plot.mp.cv.val(cv.dev, NA, ...)
  cv.dev
}

#' @author Valentijn de Jong
#' @method print   mp.cv.val
#' @export
print.mp.cv.val <- function(x, ...) {
  if (!is.null(x[["perf"]])) {
    cat("Cross-validation at stratum level yields the following performance: \n")
    print(x[["perf"]][ , 1:7])
  }
  if (!is.null(x[["gen.all"]]))
    cat("\n")
  cat("Generalizability:\n")
  print(x[["gen.all"]])
  cat("\n")
}

# Make a new meta-model to be cross-validated for metapred
# formula formula
# data data.frame
# estFUN function for estimation, e.g glm
# st.i vector, cluster indicators. char or numeric
# st.u unique values in st.i
# folds list, folds as in utils
# two.stage logical, is it a two stage method?
# estFUN function, for estimating, e.g. glm
# metaFUN function, for producing meta model
# meta.method character, option for metaFUN
# Returns mp.cv.dev, which is a list of meta-analytic models developed on dev folds
mp.cv.dev <- function(formula, data, st.i, st.u, folds, 
                      estFUN = glm, metaFUN = urma, meta.method = "DL", change = NULL, two.stage = TRUE, ...) {
  out <- list(...) # possibly contains family
  if (!is.null(out$family)) {
    if (is.character(out$family)) # Ensure that family is the output of family(). Taken directly from glm.
      out$family <- get(out$family, mode = "function", envir = parent.frame())
    if (is.function(out$family)) 
      out$family <- out$family()
    if (is.null(out$family$family)) {
      print(family)
      stop("'family' not recognized")
    }
  }
  out[["formula"]] <- formula
  out[["changed"]] <- change
  
  if (two.stage) {
    ## Cluster part 
    out[["stratified.fit"]] <- mp.stratified.fit(formula = formula, data = data, st.i = st.i, st.u = st.u, 
                                                 estFUN = estFUN, ...)  
    out[["n.clusters"]] <- length(out[["stratified.fit"]])
    
    ## cv meta-part
    out[["cv"]] <- mp.cv.meta.fit(stratified.fit = out[["stratified.fit"]], folds = folds, 
                                  metaFUN = metaFUN, meta.method = meta.method)
  } else {
    out[["cv"]] <- mp.1st.fit(formula = formula, data = data, st.i = st.i, st.u = st.u, folds = folds,
                              estFUN = estFUN, ...)
  }
  
  out[["n.cv"]] <- length(out[["cv"]])
  class(out) <- c("mp.cv.dev") 
  out
}

# For some reason I have to manually force it to find the right S3 methods.
#' @author Valentijn de Jong
#' @method print   mp.cv.dev
#' @export
print.mp.cv.dev <- function(x, ...) {
  if (!is.null(x[["stratified.fit"]]))
    print.mp.stratified.fit(x[["stratified.fit"]]) 
  cat("\n")
  print.mp.cv.meta.fit(x[["cv"]])
}

# Obtain fitted values of a mp.cv.dev object
# object mp.cv.dev object.
# Returns list of predicted values, with length equal to number of strata.
#' @author Valentijn de Jong
#' @method fitted   mp.cv.dev
#' @export
fitted.mp.cv.dev <- function(object, data, folds, st.i, predFUN = NULL, two.stage = TRUE, ...) {
  predictMethod <- getPredictMethod(fit = object, two.stage = two.stage, predFUN = predFUN, ...)
  p <- list()
  
  for (i in seq_len(object[["n.cv"]]))
    p[[i]] <- predictMethod(object = object, 
                            newdata = data[folds[["val"]][[i]] == st.i, ], 
                            type = "response",
                            b = coef(object[["cv"]][[i]]),
                            f = formula(object), 
                            two.stage = two.stage, 
                            ...)
  # If some other fold method than l1o is applied, this may otherwise cause an error.
  # better method TBI.
  if (all(unlist(lapply(folds$val, length)) == 1))
    names(p) <- unlist(folds$val)
  p
}

#' @author Valentijn de Jong
#' @method family   mp.cv.dev
#' @export
family.mp.cv.dev <- function(object, ...)
  object$family

# Estimate a one-stage or non-stratified model on the develoment (!) strata.
mp.1st.fit <- function(formula, data, st.i, st.u, folds, estFUN, ...) {
  out <- list()
  
  for (fo in seq_along(folds[["dev.i"]])) {
    fit <- estFUN(formula = formula, 
                  data = data[st.i %in% folds[["dev"]][[fo]], ],
                  ...)
    out[[getclName(folds[["dev"]][[fo]])]] <- mp.stratum.fit(fit) 
  }
  
  class(out) <- c("mp.1st.fit", class(out))
  out
}

# Recalibrate a mp.meta.fit
# mp.meta.fit mp.meta.fit
# newdata # ne data.frame, containing only the validation sample (unlike this example!)
# formula, original formula
# estFUN estimation function, e.g. glm
# Returns same object with updated coefficients.
mp.recal.meta.fit <- function(meta.fit, formula, newdata, estFUN, ...) {
  meta.fit[["formula"]] <- formula
  meta.fit[["orig coef"]] <- coef(meta.fit)
  meta.fit[["coefficients"]] <- recalibrate(meta.fit, newdata = newdata, estFUN = estFUN, ...)[["coefficients"]]
  meta.fit[["formula"]] <- NULL ### necessary ??
  meta.fit
}

# Recalibrate a mp.cv.dev
# mp.cv.dev mp.cv.dev
# newdata data.frame
# formula, original formula
# estFUN estimation function, e.g. glm
# folds folds list, as in utils.
# Returns same object with updated coefficients.
mp.cv.recal <- function(cv.dev, newdata, folds, estFUN) {
  for (i in seq_along(cv.dev[["cv"]]))
    cv.dev[["cv"]][[i]] <- mp.recal.meta.fit(meta.fit = cv.dev[["cv"]][[i]], 
                                             formula = cv.dev[["formula"]],
                                             newdata = newdata[folds[["val.i"]][[i]] == i, ], 
                                             estFUN = estFUN, 
                                             family = if (!is.null(cv.dev$family)) cv.dev$family else NULL)
  cv.dev
}

# Make new cv.models for mp.meta.fit for metapred
# This function selects the right fold for development, and calls the fitting function.
# stratified.fit list of mp.stratum.fit objects
# folds list of fold divisions, as given by l1o, or bootstrap in the utils.
# metaFUN function for estimating meta-analytic models, e.g. urma (in utils)
# meta.method options for metaFUN
# Returns an object of class mp.cv.meta.fit, which is a list of meta-analytic prediction models
mp.cv.meta.fit <- function(stratified.fit, folds, metaFUN = urma, meta.method = "DL") {
  out <- list()
  
  for (fo in seq_along(folds[["dev.i"]]))
    out[[getclName(folds[["dev"]][[fo]])]] <- mp.meta.fit(stratified.fit = stratified.fit[folds[["dev.i"]][[fo]]],
                                                          metaFUN = metaFUN, meta.method = meta.method)
  
  class(out) <- c("mp.cv.meta.fit", stratified.fit[[1]]$stratum.class)
  out
}

#' @author Valentijn de Jong
#' @method coef   mp.cv.meta.fit
#' @export
coef.mp.cv.meta.fit <- function(object, ...) 
  t(as.data.frame(lapply(object, `[[`, "coefficients"), check.names = FALSE))

#' @author Valentijn de Jong
#' @method print   mp.cv.meta.fit
#' @export
print.mp.cv.meta.fit <- function(x, ...) {
  cat("Meta-analytic models, estimated in", length(x), "fold combinations. Coefficients: \n")
  print(coef(x))
  
  if (!is.null(x[[1]][["orig coef"]])) {
    cat("\n")
    cat("Original coefficients before recalibration: \n")
    print(t(as.data.frame(lapply(x, `[[`, "orig coef"), check.names = FALSE) ))
  }
}

# Make new meta model (i.e. model fitted on multiple clusters) for ?? for metapred
# stratified.fit mp.stratified.fit
# metaFUN function for estimating meta-analytic models, e.g. urma (this file)
# meta.method options for metaFUN
# Returns object of class mp.cv.model, which is a meta-analytic prediction model
mp.meta.fit <- function(stratified.fit, metaFUN = urma, meta.method = "DL") {
  out <- list()
  
  b <- coef.mp.stratified.fit(stratified.fit) # again i have to point it to the right method.
  variances <- variances.mp.stratified.fit(stratified.fit)
  vcov <- vcov.mp.stratified.fit(stratified.fit)
  
  meta <- metaFUN(coefficients = b, variances = variances, vcov = vcov, method = meta.method) 
  
  out[["coefficients"]] <- meta$coefficients
  out[["variances"]]    <- meta$variances
  out[["se"]]           <- meta$se
  out[["ci.lb"]]        <- meta$ci.lb
  out[["ci.ub"]]        <- meta$ci.ub
  out[["tau"]]          <- meta$tau
  out[["tau2"]]         <- meta$tau2
  out[["se.tau2"]]      <- meta$se.tau2
  out[["pi"]]           <- c(meta$pi.lb, meta$pi.ub)
  out[["pi.lb"]]        <- meta$pi.lb
  out[["pi.ub"]]        <- meta$pi.ub
  out[["nobs.strata"]]  <- sapply(stratified.fit, nobs)
  out[["nobs"]]         <- sum(out[["nobs.strata"]])
  
  class(out) <- "mp.meta.fit"
  out
}

#' @author Valentijn de Jong
#' @method print   mp.meta.fit
#' @export
print.mp.meta.fit <- function(x, ...) {
  cat("Coefficients: ", "\n")
  print(coef(x))
}

# Make new stratified models for mp.meta.fit for metapred
# formula formula
# data data.frame
# st.i vector, strata indicators. char or numeric
# st.u unique values in st.i
# estFUN function for estimation, e.g glm
# ... other, e.g. family for glm.
# Returns mp.stratified.fit
mp.stratified.fit <- function(formula, data, st.i, st.u, estFUN, ...) {
  out <- list()
  
  # for (st.this in st.u)
  #   out[[getclName(st.this)]] <- mp.stratum.fit(estFUN(formula = formula, data = data[st.i == st.this, ], ...) )
  
  # Replace with: (?)
  # for (st.this in sort(unique(unlist(folds[["dev"]]) )) )
  for (st.this in st.u)
    out[[getclName(st.this)]] <- mp.stratum.fit(estFUN(formula = formula, data = data[st.i %in% st.this, ], ...) )
  
  class(out) <- "mp.stratified.fit"
  out
}

#' @author Valentijn de Jong
#' @method coef   mp.stratified.fit
#' @export
coef.mp.stratified.fit <- function(object, ...)
  as.data.frame(t(as.data.frame(lapply(object, `[[`, "coefficients", drop = FALSE))))

#' @author Valentijn de Jong
#' @method vcov   mp.stratified.fit
#' @export
vcov.mp.stratified.fit <- function(object, ...) {
  l <- lapply(object, vcov)
  array(unlist(l), dim = c(nrow(l[[1L]]), ncol(l[[1L]]), length(l)), 
        dimnames = list(rownames(l[[1L]]), colnames(l[[1L]]), names(l)))
}

#' @author Valentijn de Jong
#' @method vcov   mp.stratum.fit
#' @export
vcov.mp.stratum.fit <- function(object, ...)
  object$vcov

#' @author Valentijn de Jong
#' @method print   mp.stratified.fit
#' @export
print.mp.stratified.fit <- function(x, ...) {
  cat("Prediction models estimated in", length(x), "strata. Coefficients:\n")
  print(coef(x))
}

# Make a new mp.stratum.fit for mp.stratified.fit for metapred.
# fit model fit object, e.g. glm object
# Returns mp.stratum.fit
mp.stratum.fit <- function(fit) {
  out <- list()
  out[["coefficients"]] <- getCoefs(fit)
  out[["variances"]]    <- getVars(fit)
  # out[["covar"]]        <- getCoVars(fit)
  out[["vcov"]]         <- getCoVars(fit)
  out[["nobs"]]         <- nobs(fit, use.fallback = TRUE)
  out[["stratum.class"]]<- class(fit)
  
  class(out) <- "mp.stratum.fit"
  out
}

#' @author Valentijn de Jong
#' @method print   mp.stratum.fit
#' @export
print.mp.stratum.fit <- function(x, ...) {
  cat("Coefficients: ", "\n")
  print(coef(x))
}

# Too many args that might not be available for this method. I have abandoned it.
# predict.mp.stratum.fit <- function(object, newdata = NULL, family = NULL, 
# formula = NULL, method = NULL, type = "response") {
#   if (!is.null(family))
#     object$family <- family
#   if (is.null(formula))
#     stop("predict.mp.stratum.fit needs a formula")
#   else
#     object$formula <- formula
#   if (is.null(method))
#     stop("predict.mp.stratum.fit needs a prediction method")
# 
#   method(object = object, newdata = newdata, type = type)
# }

# As I would have had to implement it a million times:
family.default <- function(object, ...) 
  object$family

#' Standard errors and variances
#' 
#' Obtain standard errors or variances of a model fit
#' 
#' @aliases variances se tau tau2
#' 
#' @author Valentijn de Jong
#' 
#' @usage se(object, ...)
#' variances(object, ...)
#' tau(object, ...)
#' tau2(object, ...)
#' 
#' @param object A model fit object
#' @param ... other arguments
#' 
#' @return For \code{se} the standard errors of \code{object}, and for 
#' \code{variances} the variances. For \code{tau} the heterogeneity of the coefficients.
#' @export
se <- function(object, ...)
  UseMethod("se", object)

#' @export
se.default <- function(object, ...)
  sqrt(variances(object, ...))

#' @export
variances <- function(object, ...)
  UseMethod("variances", object)

#' @export
variances.default <- function(object, ...) {
  if (!is.null(v <- object$var)) # Non-exact match intented :)
    return(v)
  return(diag(vcov(object, ...))) # nrow and ncol are optional arguments! Ignore warning.
}

#' @export
variances.mp.perf <- function(object, ...)
  object$variances

#' @export
# variances.metapred <- function(object, ...)
#   variances(object[["global.model"]])

variances.metapred <- function(object, select = "global", ...)
  variances(subset(object, select = select, ...))

#' @export
variances.mp.cv <- function(object, select = "cv", ...) 
  variances(object[[select, exact = FALSE]])

#' @export
variances.mp.stratified.fit <- function(object, ...)
  t(as.data.frame(lapply(object, variances), check.names = FALSE))

#' @export
variances.mp.1st.fit <- function(object, ...) # Experimental
  t(as.data.frame(lapply(object, variances), check.names = FALSE))

#' @export
variances.mp.cv.meta.fit <- function(object, ...)
  variances.mp.stratified.fit(object, ...)

#' @export
variances.auc <- function(object, ...) 
  pROC::var(object)

#' @export
tau2 <- function(object, ...)
  tau(object, ...)^2

#' @export
tau <- function(object, ...)
  UseMethod("tau")

#' @export
tau.metapred <- function(object, ...)
  tau(subset(object, ...), ...)

#' @export
tau.default <- function(object, ...)
  object$tau

#' @author Valentijn de Jong
#' @method ci   listofperf
#' @importFrom pROC ci
#' @export
ci.listofperf <- function(object, ...) {
  # object <- list(...)$object  
  if (inherits(object[[1]], "recal")) {
    z <- lapply(object, ci.recal)
    return(data.frame(
      theta       = sapply(z, `[[`, 2),
      theta.ci.lb = sapply(z, `[[`, 1),
      theta.ci.ub = sapply(z, `[[`, 3)
    ))
  }
  if (inherits(object[[1]], "lm")) {
    z <- lapply(object, confint)
    return(data.frame(theta       = sapply(object, `[[`, 1),
                      theta.ci.lb = sapply(z, `[[`, 1),
                      theta.ci.ub = sapply(z, `[[`, 2)
    ))
  }
  if (inherits(object[[1]], "auc")) {
    z <- lapply(object, pROC::ci)
    return(data.frame(
      theta       = sapply(z, `[[`, 2),
      theta.ci.lb = sapply(z, `[[`, 1),
      theta.ci.ub = sapply(z, `[[`, 3)
    ))
  }
  if (inherits(object[[1]], "mse")) {
    z <- lapply(object, ci.mse)
    return(data.frame(
      theta       = sapply(z, `[[`, 2),
      theta.ci.lb = sapply(z, `[[`, 1),
      theta.ci.ub = sapply(z, `[[`, 3)
    ))
  }
  stop(paste("ci.listofperf does not recognize objects of class", class(object)))
}
# confint.listofperf <- function(object, ...) {
#   if (inherits(object[[1]], "lm")) {
#     z <- lapply(object, confint, parm = parm, level = level)
#     return(data.frame(theta       = sapply(object, `[[`, 1),
#                       theta.ci.lb = sapply(z, `[[`, 1),
#                       theta.ci.ub = sapply(z, `[[`, 2)
#                       ))
#   }
#   if (inherits(object[[1]], "auc")) {
#     z <- lapply(object, pROC::ci)
#     return(data.frame(
#       theta       = sapply(z, `[[`, 2),
#       theta.ci.lb = sapply(z, `[[`, 1),
#       theta.ci.ub = sapply(z, `[[`, 3)
#     ))
#   }
#   stop("confint.listofperf does not recognize the object.")
# }

# #' @importFrom pROC ci
# #' @export
# confint.auc <- function(object, parm, level, ...) {
#   
# }

#' @export
ci.recal <- function(object, conf = .95, ...) { # To be implemented in some other class. Temporarily for lm.
  ses <- se(object, ...)
  coefs <- coef(object, ...)
  z <- qt(1 - (1 - conf)/2, df = object$df.residual) # z = t distributed
  data.frame("ci.lb" = coefs - z * ses, "estimate" = coefs, "ci.ub" = coefs + z * ses)
}

#' @export
ci.mse <- function(object, conf = .95, ...) {
  ses <- se(object, ...)
  est <- object[["estimate"]]
  z <- qt(1 - (1 - conf)/2, df = object$n - 1) # z = t distributed
  data.frame("ci.lb" = est - z * ses, "estimate" = est, "ci.ub" = est + z * ses)
}


#' Generalizability estimates
#' 
#' Obtain generalizability estimates from a model fit.
#' 
#' @aliases gen generalizability
#' 
#' @author Valentijn de Jong
#' 
#' @usage gen(object, ...)
#' generalizability(object, ...)
#' 
#' @param object A model fit object, either a \link{metapred} or \code{subset(metapred)} object.
#' @param ... By default, the final model is selected. This parameter allows other arguments to be passed 
#' to \link{subset.metapred} such that the generalizability estimates of other steps/models may be
#' returned.. 
#' 
#' @export
gen <- function(object, ...) 
  UseMethod("gen")

#' @export
generalizability <- gen

#' @export
gen.metapred <- function(object, genFUN = 1, ...) 
  gen(subset(object, ...), genFUN = genFUN, ...)

#' @export
gen.mp.cv.val <- function(object, genFUN = 1, ...)
  object$gen.all[[genFUN]]


#' Performance estimates
#' 
#' Obtain performance estimates from a model fit.
#' 
#' @aliases perf performance 
#' 
#' @author Valentijn de Jong
#' 
#' @usage perf(object, ...)
#' performance(object, ...)
#' 
#' @param object A model fit object, either a \link{metapred} or \code{subset(metapred)} object.
#' @param ... By default, the final model is selected. This parameter allows other arguments to be
#' passed to \link{subset.metapred} such that the performance estimates of other steps/models may be
#' returned.. 
#' 
#' @export
perf <- function(object, ...) 
  UseMethod("perf")

#' @export
performance <- perf

#' @export
perf.metapred <- function(object, perfFUN = 1, ...) 
  perf(subset(object, ...), perfFUN = perfFUN, ...)

#' @export
perf.mp.cv.val <- function(object, perfFUN = 1, ...)
  object[["perf.all"]][[perfFUN]]