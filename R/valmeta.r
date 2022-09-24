#' Meta-analysis of prediction model performance
#'
#' This function provides summary estimates for the concordance statistic, 
#' the total observed-expected ratio or the calibration slope. Where 
#' appropriate, data transformations are applied and missing information 
#' is derived from available quantities. Unless specified otherwise, all 
#' meta-analysis models assume random effects and are fitted using restricted 
#' maximum likelihood estimation with the \pkg{metafor} package 
#' (Viechtbauer 2010). Further, confidence intervals for the average 
#' performance are based on the Hartung-Knapp-Sidik-Jonkman method. 
#' When conducting a Bayesian meta-analysis, the R packages \pkg{runjags} and 
#' \pkg{rjags} must be installed.
#' 
#' @param  measure A character string indicating which summary performance 
#' measure should be calculated. Options are
#' \code{"cstat"} (meta-analysis of the concordance statistic) and \code{"OE"} 
#' (meta-analysis of the total observed-expected ratio). See `Details' for more information.
#' @param cstat Optional vector with the estimated c-statistic for each valiation
#' @param cstat.se Optional vector with the standard error of the estimated c-statistics
#' @param cstat.cilb Optional vector to specify the lower limits of the confidence interval.
#' @param cstat.ciub Optional vector to specify the upper limits of the confidence interval.
#' @param cstat.cilv Optional vector to specify the levels of aformentioned confidence interval limits. 
#' (default: 0.95, which corresponds to the 95\% confidence interval).
#' @param sd.LP Optional vector with the standard deviation of the linear predictor (prognostic index)
#' @param OE Optional vector with the estimated ratio of total observed versus total expected events
#' @param OE.se Optional vector with the standard errors of the estimated O:E ratios
#' @param OE.cilb Optional vector to specify the lower limits of the confidence interval for \code{OE}.
#' @param OE.ciub Optional vector to specify the upper limits of the confidence interval for \code{OE}.
#' @param OE.cilv Optional vector to specify the levels of aformentioned confidence interval limits. 
#' (default: 0.95, which corresponds to the 95\% confidence interval).
#' @param citl Optional vector with the estimated calibration-in-the-large for each valiation
#' @param citl.se Optional vector with the standard error of the estimated calibration-in-the-large statistics
#' @param N Optional vector with the total number of participants for each valiation
#' @param O Optional vector with the total number of observed events for each valiation
#' (if specified, during time \code{t.val})
#' @param E Optional vector with the total number of expected events for each valiation 
#' (if specified, during time \code{t.val})
#' @param Po Optional vector with the (cumulative) observed event probability for each valiation
#' (if specified, during time \code{t.val})
#' @param Po.se Optional vector with the standard errors of \code{Po}.
#' @param Pe Optional vector with the (cumulative) expected event probability for each validation
#' (if specified, during time \code{t.val})
#' @param data optional data frame containing the variables given to the arguments above.
#' @param method Character string specifying whether a fixed- or a random-effects model should be fitted. 
#' A fixed-effects model is fitted when using \code{method="FE"}. Random-effects models are fitted by setting method 
#' equal to one of the following: \code{"REML"} (Default), \code{"DL"}, \code{"HE"}, \code{"SJ"}, \code{"ML"}, 
#' \code{"EB"}, \code{"HS"}, \code{"GENQ"} or \code{"BAYES"}. See 'Details'.
#' @param test Optional character string specifying how test statistics and confidence intervals for the fixed effects 
#' should be computed. By default (\code{test="knha"}), the method by Knapp and Hartung (2003) is used for 
#' adjusting test statistics and confidence intervals. Type '\code{?rma}' for more details.
#' @param verbose If TRUE then messages generated during the fitting process will be displayed.
#' @param slab Optional vector specifying the label for each study
#' @param n.chains Optional numeric specifying the number of chains to use in the Gibbs sampler 
#' (if \code{method="BAYES"}). More chains will improve the sensitivity of the convergence diagnostic, but will 
#' cause the simulation to run more slowly. The default number of chains is 4.
#' @param pars A list with additional arguments.  See 'Details' for more information. The following parameters configure the MCMC sampling procedure:  
#' \code{hp.mu.mean} (mean of the prior distribution of the random effects model, defaults to 0), 
#' \code{hp.mu.var} (variance of the prior distribution of the random effects model, defaults to 1000), 
#' \code{hp.tau.min} (minimum value for the between-study standard deviation, defaults to 0), 
#' \code{hp.tau.max} (maximum value for the between-study standard deviation, defaults to 2), 
#' \code{hp.tau.sigma} (standard deviation of the prior distribution for the between-study standard-deviation), 
#' \code{hp.tau.dist} (prior distribution for the between-study standard-deviation. Defaults to \code{"dunif"}), 
#' \code{hp.tau.df} (degrees of freedom for the prior distribution for the between-study standard-deviation. 
#' Defaults to 3). Other arguments are \code{method.restore.c.se} (method for restoring missing estimates for the standard error 
#' of the c-statistic. See \code{\link{ccalc}} for more information), 
#' \code{model.cstat} (The likelihood/link for modeling the c-statistic; see "Details"), 
#' \code{model.oe} (The likelihood/link for modeling the O:E ratio; see "Details"),
#' \code{seed} (an integer to indicate the random seed).
#' @param \ldots Additional arguments that are passed to \pkg{rma} or \pkg{runjags} (if \code{method="BAYES"}).
#' 
#' @details 
#' \subsection{Meta-analysis of the concordance statistic}{
#' A summary estimate for the concorcance (c-) statistic can be obtained by specifying \code{measure="cstat"}.
#' The c-statistic is a measure of discrimination, and indicates the ability of a prediction model to 
#' distinguish between patients developing and not developing the outcome. The c-statistic typically ranges 
#' from 0.5 (no discriminative ability) to 1 (perfect discriminative ability). When missing, the c-statistic 
#' and/or its standard error are derived from other reported information. 
#' See \code{\link{ccalc}} for more information.
#' 
#' By default, it is assumed that the logit of the c-statistic is Normally distributed within and across studies
#' (\code{pars$model.cstat = "normal/logit"}). Alternatively, it is possible to assume that the raw c-statistic 
#' is Normally distributed across studies \code{pars$model.cstat = "normal/identity"}. 
#' } 
#' 
#' \subsection{Meta-analysis of the total observed versus expected ratio}{
#' A summary estimate for the total observed versus expected (O:E) ratio can be 
#' obtained by specifying \code{measure="OE"}. The total O:E ratio provides a 
#' rough indication of the overall model calibration (across the entire range 
#' of predicted risks). When missing, the total O:E ratio and/or its standard 
#' error are derived from other reported information. See \code{\link{oecalc}} 
#' for more information.
#' 
#' For frequentist meta-analysis, within-study variation can either be modeled 
#' using a Normal (\code{model.oe = "normal/log"} or 
#' \code{model.oe = "normal/identity"}) or a Poisson distribution 
#' (\code{model.oe = "poisson/log"}). 
#' 
#' When performing a Bayesian meta-analysis, all data are modeled using a 
#' one-stage random effects (hierarchical related regression) model.
#' In particular, a binomial distribution (if \code{O}, \code{E} and 
#' \code{N} is known), a Poisson distribution (if only \code{O} and 
#' \code{E} are known) or a Normal distribution (if \code{OE} and 
#' \code{OE.se} or \code{OE.95CI} are known) is selected separately for 
#' each study.
#' 
#' }
#' 
#' \subsection{Bayesian meta-analysis}{
#' All Bayesian meta-analysis models assume the presence of random effects. 
#' Summary estimates are based on the posterior mean. Credibility and 
#' prediction intervals are directly obtained from the corresponding posterior 
#' quantiles.
#' 
#' The prior distribution for the (transformed) performance estimate is modeled 
#' using a Normal distribution, with mean \code{hp.mu.mean} (defaults to 0) 
#' and variance \code{hp.mu.var} (defaults to 1000). 
#' For meta-analysis of the total O:E ratio, the maximum value for 
#' \code{hp.mu.var} is 100.
#' 
#' By default, the prior distribution for the between-study standard deviation 
#' is modeled using a uniform distribution (\code{hp.tau.dist="dunif"}), 
#' with boundaries \code{hp.tau.min} and \code{hp.tau.max}. Alternative choices 
#' are a truncated Student-t distribution (\code{hp.tau.dist="dhalft"}) with a 
#' mean of \code{hp.tau.mean}, a standard deviation of \code{hp.tau.sigma} and 
#' \code{hp.tau.df} degrees of freedom. This distribution is again restricted 
#' to the range \code{hp.tau.min} to \code{hp.tau.max}.
#' }
#' 
#' @note The width of calculated confidence, credibility and prediction 
#' intervals can be specified using \code{level} in the \code{pars} argument 
#' (defaults to 0.95).
#' 
#' @return An object of class \code{valmeta} with the following elements:
#' \describe{
##'  \item{"data"}{array with (transformed) data used for meta-analysis, 
##'  and method(s) used for restoring missing information. }
##'  \item{"measure"}{character string specifying the performance measure that 
##'  has been meta-analysed.}
##'  \item{"method"}{character string specifying the meta-analysis method.}
##'  \item{"model"}{character string specifying the meta-analysis model 
##'  (link function).}
##'  \item{"est"}{summary estimate for the performance statistic. For Bayesian 
##'  meta-analysis, the posterior median is returned.}
##'  \item{"ci.lb"}{lower bound of the confidence (or credibility) interval of 
##'  the summary performance estimate.}
##'  \item{"ci.ub"}{upper bound of the confidence (or credibility) interval of 
##'  the summary performance estimate.}
##'  \item{"pi.lb"}{lower bound of the (approximate) prediction interval of the 
##'  summary performance estimate.}
##'  \item{"pi.ub"}{upper bound of the (approximate) prediction interval of the 
##'  summary performance estimate.}
##'  \item{"fit"}{the full results from the fitted model. }
##'  \item{"slab"}{vector specifying the label of each study.}
##' }
#' @references 
#' Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. 
#' A guide to systematic review and meta-analysis of prediction model 
#' performance. \emph{BMJ}. 2017;\bold{356}:i6460. \doi{10.1136/bmj.i6460}
#' 
#' Debray TPA, Damen JAAG, Riley R, Snell KIE, Reitsma JB, Hooft L, et al. 
#' A framework for meta-analysis of prediction model studies with binary and 
#' time-to-event outcomes. \emph{Stat Methods Med Res}. 2019;\bold{28}:2768--86. 
#' \doi{10.1177/0962280218785504}
#' 
#' Riley RD, Tierney JF, Stewart LA. Individual participant data meta-analysis: 
#' a handbook for healthcare research. Hoboken, NJ: Wiley; 2021. 
#' ISBN: 978-1-119-33372-2.
#' 
#' Steyerberg EW, Nieboer D, Debray TPA, van Houwelingen HC. Assessment of 
#' heterogeneity in an individual participant data meta-analysis of prediction 
#' models: An overview and illustration. \emph{Stat Med}. 2019;
#' \bold{38}:4290--309. \doi{10.1002/sim.8296}
#' 
#' Viechtbauer W. Conducting Meta-Analyses in R with the metafor Package. 
#' \emph{Journal of Statistical Software}. 2010; \bold{36}. 
#' \doi{10.18637/jss.v036.i03}
#'   
#' @seealso \code{\link{ccalc}} to calculate concordance statistics 
#' and corresponding standard errors, \code{\link{oecalc}} to calculate the 
#' total O:E ratio and corresponding standard errors, 
#' \code{\link{plot.valmeta}} to generate forest plots
#' 
#' @examples 
#' ######### Validation of prediction models with a binary outcome #########
#' data(EuroSCORE)
#' 
#' # Meta-analysis of the c-statistic (random effects)
#' fit <- valmeta(cstat=c.index, cstat.se=se.c.index, cstat.cilb=c.index.95CIl, 
#'                cstat.ciub=c.index.95CIu, cstat.cilv=0.95, N=n, O=n.events, 
#'                slab=Study, data=EuroSCORE)
#' plot(fit)
#' 
#' # Nearly identical results when we need to estimate the SE
#' valmeta(cstat=c.index,  N=n, O=n.events, slab=Study, data=EuroSCORE)
#' 
#' # Two-stage meta-analysis of the total O:E ratio (random effects)
#' valmeta(measure="OE", O=n.events, E=e.events, N=n, slab=Study, data=EuroSCORE)    
#' valmeta(measure="OE", O=n.events, E=e.events, data=EuroSCORE)       
#' valmeta(measure="OE", Po=Po, Pe=Pe, N=n, data=EuroSCORE)
#' 
#' \dontrun{
#' # One-stage meta-analysis of the total O:E ratio (random effects)
#' valmeta(measure="OE", O=n.events, E=e.events, data=EuroSCORE, method="ML", 
#'         pars=list(model.oe="poisson/log"))
#' 
#' # Bayesian random effects meta-analysis of the c-statistic
#' fit2 <- valmeta(cstat=c.index, cstat.se=se.c.index, cstat.cilb=c.index.95CIl,
#'                 cstat.ciub=c.index.95CIu, cstat.cilv=0.95, N=n, O=n.events, 
#'                 data=EuroSCORE, method="BAYES", slab=Study)
#' 
#' # Bayesian one-stage random effects meta-analysis of the total O:E ratio
#' # Consider that some (but not all) studies do not provide information on N
#' # A Poisson distribution will be used for studies 1, 2, 5, 10 and 20
#' # A Binomial distribution will be used for the remaining studies
#' EuroSCORE.new <- EuroSCORE
#' EuroSCORE.new$n[c(1, 2, 5, 10, 20)] <-  NA
#' pars <- list(hp.tau.dist="dhalft",   # Prior for the between-study standard deviation
#'              hp.tau.sigma=1.5,       # Standard deviation for 'hp.tau.dist'
#'              hp.tau.df=3,            # Degrees of freedom for 'hp.tau.dist'
#'              hp.tau.max=10,          # Maximum value for the between-study standard deviation
#'              seed=5)                 # Set random seed for the simulations
#' fit3 <- valmeta(measure="OE", O=n.events, E=e.events, N=n, data=EuroSCORE.new,
#'         method="BAYES", slab=Study, pars=pars)
#' plot(fit3)
#' print(fit3$fit$model) # Inspect the JAGS model
#' print(fit3$fit$data)  # Inspect the JAGS data
#' } 
#' 
#' ######### Validation of prediction models with a time-to-event outcome #########
#' data(Framingham)
#' 
#' # Meta-analysis of total O:E ratio after 10 years of follow-up
#' valmeta(measure="OE", Po=Po, Pe=Pe, N=n, data=Framingham)
#' 
#' @keywords meta-analysis discrimination  calibration
#' 
#' 
#' @export
#' @import metafor
#' @import mvtnorm
#' @importFrom dplyr select mutate %>%
#' @importFrom lme4 glmer
#' @importFrom stats coef coefficients dnorm glm nobs optim pchisq qnorm qt pt rnorm runif confint poisson
#' predict vcov as.formula formula model.frame model.frame.default update.formula family
valmeta <- function(measure="cstat", cstat, cstat.se, cstat.cilb, cstat.ciub, cstat.cilv,
                    sd.LP, OE, OE.se, OE.cilb, OE.ciub, OE.cilv, citl, citl.se, N, O, E, Po, Po.se, Pe, data, 
                    method="REML", test="knha", verbose=FALSE, slab, n.chains = 4, pars, ...) {
  
  pars.default <- .initiateDefaultPars(pars, type = "valmeta")
  
  ### check if data argument has been specified
  if (missing(data))
    data <- NULL
  
  ### need this at the end to check if append=TRUE can actually be done
  no.data <- is.null(data)
  
  ### check if data argument has been specified
  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data))
      data <- data.frame(data)
  }
  
  #######################################################################################
  # Retrieve all data
  #######################################################################################
  mf <- match.call()
  
  mf.slab       <- mf[[match("slab",   names(mf))]]
  slab          <- eval(mf.slab,   data, enclos = sys.frame(sys.parent()))
  mf.cstat      <- mf[[match("cstat", names(mf))]]
  cstat         <- eval(mf.cstat, data, enclos = sys.frame(sys.parent()))
  mf.cstat.se   <- mf[[match("cstat.se", names(mf))]]
  cstat.se      <- eval(mf.cstat.se, data, enclos = sys.frame(sys.parent()))
  mf.cstat.cilb <- mf[[match("cstat.cilb", names(mf))]]
  cstat.cilb    <- eval(mf.cstat.cilb, data, enclos = sys.frame(sys.parent()))
  mf.cstat.ciub <- mf[[match("cstat.ciub", names(mf))]]
  cstat.ciub    <- eval(mf.cstat.ciub, data, enclos = sys.frame(sys.parent()))
  mf.cstat.cilv <- mf[[match("cstat.cilv", names(mf))]]
  cstat.cilv    <- eval(mf.cstat.cilv, data, enclos = sys.frame(sys.parent()))
  mf.sd.LP      <- mf[[match("sd.LP", names(mf))]]
  sd.LP         <- eval(mf.sd.LP, data, enclos = sys.frame(sys.parent()))
  mf.OE         <- mf[[match("OE", names(mf))]]
  OE            <- eval(mf.OE, data, enclos = sys.frame(sys.parent()))
  mf.OE.se      <- mf[[match("OE.se", names(mf))]]
  OE.se         <- eval(mf.OE.se, data, enclos = sys.frame(sys.parent()))
  mf.OE.cilb    <- mf[[match("OE.cilb", names(mf))]]
  OE.cilb       <- eval(mf.OE.cilb, data, enclos = sys.frame(sys.parent()))
  mf.OE.ciub    <- mf[[match("OE.ciub", names(mf))]]
  OE.ciub       <- eval(mf.OE.ciub, data, enclos = sys.frame(sys.parent()))
  mf.OE.cilv    <- mf[[match("OE.cilv", names(mf))]]
  OE.cilv       <- eval(mf.OE.cilv, data, enclos = sys.frame(sys.parent()))
  mf.citl       <- mf[[match("citl", names(mf))]]
  citl          <- eval(mf.citl, data, enclos = sys.frame(sys.parent()))
  mf.citl.se    <- mf[[match("citl.se", names(mf))]]
  citl.se       <- eval(mf.citl.se, data, enclos = sys.frame(sys.parent()))
  mf.N          <- mf[[match("N", names(mf))]]
  N             <- eval(mf.N, data, enclos = sys.frame(sys.parent()))
  mf.O          <- mf[[match("O", names(mf))]]
  O             <- eval(mf.O, data, enclos = sys.frame(sys.parent()))
  mf.E          <- mf[[match("E", names(mf))]]
  E             <- eval(mf.E, data, enclos = sys.frame(sys.parent()))
  mf.Po         <- mf[[match("Po", names(mf))]]
  Po            <- eval(mf.Po, data, enclos = sys.frame(sys.parent()))
  mf.Po.se      <- mf[[match("Po.se", names(mf))]]
  Po.se         <- eval(mf.Po.se, data, enclos = sys.frame(sys.parent()))
  mf.Pe         <- mf[[match("Pe", names(mf))]]
  Pe            <- eval(mf.Pe, data, enclos = sys.frame(sys.parent()))
  
  if (is.null(method)) {
    stop("No estimation method specified!")
  }
  
  if (method == "FE") {
    test <- "z" #Do not use SJHK adjustment in a fixed effect MA
  }
  
  
  #######################################################################################
  # Check if we need to load runjags
  #######################################################################################
  if (method == "BAYES") {
    if (!requireNamespace("runjags", quietly = TRUE)) {
      stop("The package 'runjags' is currently not installed!")
    } 
    if (!requireNamespace("rjags", quietly = TRUE)) {
      stop("The package 'rjags' is currently not installed!")
    } 
    if (n.chains < 1 | n.chains %% 1 != 0) {
      stop("Invalid number of chains specified for the Gibbs sampler!")
    }
  }
  
  #######################################################################################
  # Count number of studies
  #######################################################################################
  k <- 0
  
  if (!no.data) {
    k <- dim(data)[1]
  } else if (measure == "cstat") {
    if (!is.null(cstat)) {
      k <- length(cstat)
    } else if (!is.null(cstat.se)) {
      k <- length(cstat.se)
    } else if (!is.null(cstat.cilb)) {
      k <- length(cstat.cilb)
    } else if (!is.null(cstat.ciub)) {
      k <- length(cstat.ciub)
    } else if (!is.null(sd.LP)) {
      k <- length(sd.LP)
    }
  } else if (measure == "OE") {
    if (!is.null(OE)) {
      k <- length(OE)
    } else if (!is.null(OE.se)) {
      k <- length(OE.se)
    } else if (!is.null(OE.cilb)) {
      k <- length(OE.cilb)
    } else if (!is.null(OE.ciub)) {
      k <- length(OE.ciub)
    } else if (!is.null(E)) {
      k <- length(E)
    } else if (!is.null(O)) {
      k <- length(O)
    } else if (!is.null(Po)) {
      k <- length(Po)
    } else if (!is.null(Pe)) {
      k <- length(Pe)
    } else if (!is.null(citl)) {
      k <- length(citl)
    }
  } else {
    stop("Unknown 'measure' specified.")
  }
  
  if (k < 1) stop("No data provided!")
  
  
  #######################################################################################
  # Prepare data
  #######################################################################################
  if (is.null(slab) & !no.data) {
    slab <- rownames(data)
  } else if (!is.null(slab)) {
    slab <- make.unique(as.character(slab))
  }
  
  if (is.null(cstat)) cstat <- rep(NA, times = k)
  if (is.null(cstat.se)) cstat.se <- rep(NA, times = k)
  if (is.null(cstat.cilb)) cstat.cilb <- rep(NA, times = k)
  if (is.null(cstat.ciub)) cstat.ciub <- rep(NA, times = k)
  if (is.null(cstat.cilv)) cstat.cilv <- rep(0.95, times = k)
  if (is.null(O)) O <- rep(NA, times = k)
  if (is.null(Po)) Po <- rep(NA, times = k)
  if (is.null(N)) N <- rep(NA, times = k)
  if (is.null(sd.LP)) sd.LP <- rep(NA, times = k)
  if (is.null(OE)) OE <- rep(NA, times = k)
  if (is.null(OE.se)) OE.se <- rep(NA, times = k)
  if (is.null(OE.cilb)) OE.cilb <- rep(NA, times = k)
  if (is.null(OE.ciub)) OE.ciub <- rep(NA, times = k)
  if (is.null(OE.cilv)) OE.cilv <- rep(0.95, times = k) # Assume 95% CI by default 
  if (is.null(E)) E <- rep(NA, times = k)
  if (is.null(Po.se)) Po.se <- rep(NA, times = k)
  if (is.null(Pe)) Pe <- rep(NA, times = k)
  if (missing(citl)) citl <- rep(NA, times = k)
  if (is.null(citl.se)) citl.se <- rep(NA, times = k)
  
  #######################################################################################
  # Prepare results
  #######################################################################################
  ci.lb <- ci.ub <- NA
  
  #######################################################################################
  # Prepare object
  #######################################################################################
  out <- list()
  out$call <- match.call()
  out$measure <- measure
  out$method <- method
  out$level <- pars.default$level 
  class(out) <- "valmeta"
  
  
  #######################################################################################
  # Meta-analysis of the c-statistic
  #######################################################################################
  if (measure == "cstat") {
    if (verbose) message("Extracting/computing estimates of the c-statistic ...")
    
    out$model <- pars.default$model.cstat
    
    if (out$model == "normal/identity") {
      g <- NULL
    } else if (out$model == "normal/logit") {
      g <- "log(cstat/(1-cstat))"
    } else {
      stop(paste("Meta-analysis model currently not supported: '", out$model, '"', sep = ""))
    }
    
    ds <- ccalc(cstat = cstat, 
                cstat.se = cstat.se, 
                cstat.cilb = cstat.cilb,
                cstat.ciub = cstat.ciub,
                cstat.cilv = cstat.cilv,
                sd.LP = sd.LP, 
                N = N, 
                O = O, 
                Po = Po, 
                slab = slab, 
                g = g, 
                level = pars.default$level, 
                approx.se.method = pars.default$method.restore.c.se) 
    
    ## Assign study labels
    out$slab <- rownames(ds)
    
    # Identify which studies can be used for meta-analysis
    selstudies <- which(!is.na(ds$theta) & !is.na(ds$theta.se))
    
    if (method != "BAYES") { # Use of rma
      
      # Apply the meta-analysis
      fit <- metafor::rma(yi = ds$theta, 
                          sei = ds$theta.se, 
                          method = method, 
                          test = test, 
                          slab = out$slab, ...)
      
      preds <- predict(fit, level = pars.default$level)
      
      # The predict function from metafor uses a Normal distribution for prediction intervals, 
      # Here, we will use a Student T distribution instead
      predint <- calcPredInt(coefficients(fit), 
                             sigma2 = fit$se**2, 
                             tau2 = fit$tau2, 
                             k = fit$k,
                             level = pars.default$level)
      pi.lb <- predint$lower
      pi.ub <- predint$upper
      
      ds[selstudies, "theta.blup"] <- metafor::blup(fit)$pred
      
      if (pars.default$model.cstat == "normal/logit") {
        out$est <- as.numeric(inv.logit(coefficients(fit)))
        out$ci.lb <- inv.logit(preds$ci.lb)
        out$ci.ub <- inv.logit(preds$ci.ub)
        out$pi.lb <- ifelse(method == "FE", inv.logit(ci.lb), inv.logit(pi.lb))
        out$pi.ub <- ifelse(method == "FE", inv.logit(ci.ub), inv.logit(pi.ub))
      } else if (pars.default$model.cstat == "normal/identity") {
        out$est <- as.numeric(coefficients(fit))
        out$ci.lb <- preds$ci.lb
        out$ci.ub <- preds$ci.ub
        out$pi.lb <- ifelse(method == "FE", ci.lb, pi.lb)
        out$pi.ub <- ifelse(method == "FE", ci.ub, pi.ub)
      } else {
        stop("There is no implementation for the specified meta-analysis model!")
      }
      
      out$fit <- fit
      out$numstudies <- fit$k
      out$data <- ds
    } else {
      return(run_Bayesian_REMA(call = match.call(),
                               measure = measure,
                               method = method,
                               data = ds, 
                               pars = pars.default, 
                               FUN_generate_bugs = .generateBugsCstat,
                               n.chains = n.chains, 
                               verbose = verbose, ...))
    }
    
    return(out)
  }
  #######################################################################################
  # Meta-analysis of the total OE ratio
  #######################################################################################
  if (measure == "OE") {
    if (verbose) message("Extracting/computing estimates of the total O:E ratio ...")
    

    if (pars.default$model.oe == "normal/identity") {
      g <- NULL
    } else if (pars.default$model.oe == "normal/log" |  pars.default$model.oe == "poisson/log") {
      g <- "log(OE)"
    } else {
      stop(paste("Meta-analysis model currently not supported: '", pars.default$model.oe, '"', sep = ""))
    }
    
    ds <- oecalc(OE = OE, 
                 OE.se = OE.se,
                 OE.cilb = OE.cilb,
                 OE.ciub = OE.ciub,
                 OE.cilv = OE.cilv,
                 citl = citl, 
                 citl.se = citl.se,
                 N = N,
                 O = O,
                 E = E, 
                 Po = Po,
                 Po.se = Po.se,
                 Pe = Pe, 
                 slab = slab, 
                 g = g,
                 level = pars.default$level) 
    
    
    if (method == "BAYES") {
      if (verbose) print("Performing Bayesian one-stage meta-analysis...")
      return(run_Bayesian_MA_oe(call = match.call(),
                                measure = measure,
                                method = method,
                                data = ds, 
                                pars = pars.default, 
                                n.chains = n.chains, 
                                verbose = verbose, ...))

    }
    
    if (method != "BAYES") { # Use of rma
      ## Assign study labels
      out$slab <- rownames(ds)
      out$model <- pars.default$model.oe
      
      if (pars.default$model.oe == "normal/identity") {
        if (verbose) print("Performing a frequentist two-stage meta-analysis...")
        fit <- metafor::rma(yi = ds$theta, sei = ds$theta.se, data = ds, method = method, test = test, slab = out$slab, ...) 
        preds <- predict(fit, level = pars.default$level)
        
        # The predict function from metafor uses a Normal distribution for prediction intervals, 
        # Here, we will use a Student T distribution instead
        predint <- calcPredInt(coefficients(fit), 
                               sigma2 = fit$se**2, 
                               tau2 = fit$tau2, 
                               k = fit$k, 
                               level = pars.default$level)
        
        out$est   <- as.numeric(coefficients(fit))
        out$ci.lb <- as.numeric(preds$ci.lb)
        out$ci.ub <- as.numeric(preds$ci.ub)
        out$pi.lb <- ifelse(method == "FE", ci.lb, predint$lower)
        out$pi.ub <- ifelse(method == "FE", ci.ub, predint$upper)
        out$fit <- fit
        out$numstudies <- fit$k
        
      } else if (pars.default$model.oe == "normal/log") {
        if (verbose) print("Performing two-stage meta-analysis...")
        fit <- metafor::rma(yi = ds$theta, sei = ds$theta.se, data = ds, method = method, test = test, slab = out$slab, ...) 
        preds <- predict(fit, level = pars.default$level)
        
        # The predict function from metafor uses a Normal distribution for prediction intervals, 
        # Here, we will use a Student T distribution instead
        predint <- calcPredInt(coefficients(fit), sigma2=fit$se**2, tau2=fit$tau2, k=fit$k, level=pars.default$level)
        
        out$est    <- as.numeric(exp(coefficients(fit)))
        out$ci.lb  <- exp(preds$ci.lb)
        out$ci.ub  <- exp(preds$ci.ub)
        out$pi.lb  <- ifelse(method=="FE", exp(ci.lb), exp(predint$lower))
        out$pi.ub  <- ifelse(method=="FE", exp(ci.lb), exp(predint$upper))
        
        out$fit <- fit
        
        out$numstudies <- fit$k
      } else if (pars.default$model.oe == "poisson/log") { 
        if (verbose) print("Performing one-stage meta-analysis...")
        if (method == "ML") { 
          if (test == "knha") warning("The Sidik-Jonkman-Hartung-Knapp correction cannot be applied")
          ds$Study <- out$slab
          ds$O <- round(ds$O)
          fit <- lme4::glmer(O~1|Study, offset = log(E), family = poisson(link="log"), data = ds)
          
          # Extract the random effects
          theta.ranef <- lme4::ranef(fit, drop = TRUE, condVar = TRUE)
          #ds <- select(ds, -c("Study", "theta", "theta.se","theta.cilb", "theta.ciub"))
          
          ia <- which(!is.na(ds$O) & !is.na(ds$E))
          ds$theta.blup <- ds$theta.se.blup <- NA
          ds$theta.blup[ia] <- as.numeric(theta.ranef$Study)
          ds$theta.se.blup[ia] <- sqrt(attr(theta.ranef[[1]], "postVar"))
          
          preds.ci <- confint(fit, 
                              level = pars.default$level, 
                              quiet = !verbose, ...)
          predint <- calcPredInt(lme4::fixef(fit), 
                                 sigma2 = vcov(fit)[1,1], 
                                 tau2 = (as.data.frame(lme4::VarCorr(fit))["vcov"])[1,1], 
                                 k = lme4::ngrps(fit), 
                                 level = pars.default$level)
          
          out$est    <- as.numeric(exp(lme4::fixef(fit)))
          out$ci.lb  <- exp(preds.ci["(Intercept)",1])
          out$ci.ub  <- exp(preds.ci["(Intercept)",2])
          out$pi.lb  <- exp(predint$lower)
          out$pi.ub  <- exp(predint$upper)
          
          out$fit <- fit
          out$numstudies <- nobs(fit)
        } else if (method == "FE") { #one-stage fixed-effects meta-analysis
          fit <- glm(O~1, offset = log(E), family = poisson(link = "log"), data = ds)
          
          preds.ci <- confint(fit, level = pars.default$level, quiet = !verbose, ...)
          
          out$est   <- as.numeric(exp(coefficients(fit)))
          out$ci.lb <- out$pi.lb <- exp(preds.ci[1])
          out$ci.ub <- out$pi.ub <- exp(preds.ci[2])
          
          out$fit <- fit
          out$numstudies <- nobs(fit)
        } else {
          stop(paste("No implementation found for ", method, " estimation (model '", pars.default$model.oe, "')!", sep = ""))
        }
      } else {
        stop("Model not implemented yet!")
      }
    } 
    
    out$data <- ds
    
    return(out)
  }
}


#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method print valmeta
#' @export
print.valmeta <- function(x, ...) {
  text.stat <- attr(x$data,'estimand')
  text.model <- if (x$method == "FE") "Fixed" else "Random"
  text.ci <- if (x$method == "BAYES") "credibility" else "confidence"
  text.pi <- if (x$method == "BAYES") "" else "(approximate)"
  
  
  if (x$method != "FE") {
    cat(paste("Summary ", text.stat, " with ", x$level*100, "% ", text.ci, " and ", text.pi, " ", x$level*100, "% prediction interval:\n\n", sep=""))
    results <- c(Estimate=x$est, CIl=x$ci.lb, CIu=x$ci.ub, PIl=x$pi.lb, PIu=x$pi.ub)
  } else {
    cat(paste("Summary ", text.stat, " with ", x$level*100, "% ", text.ci, " interval:\n\n", sep = ""))
    results <- c(Estimate = x$est, CIl = x$ci.lb, CIu = x$ci.ub)
  }
  print(results)
  cat("\n")
  cat(paste("Number of studies included: ", x$numstudies))
}

#' Forest Plots
#' 
#' Function to create forest plots for objects of class \code{"valmeta"}.
#' 
#' @param x An object of class \code{"valmeta"}
#' @param \ldots Additional arguments which are passed to \link{forest}.
#' 
#' @details The forest plot shows the performance estimates of each validation with corresponding confidence 
#' intervals. A polygon is added to the bottom of the forest plot, showing the summary estimate based on the model. 
#' A 95\% prediction interval is added by default for random-effects models,  the dotted line indicates its (approximate) bounds.
#' 
#' @references 
#' Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review 
#' and meta-analysis of prediction model performance. \emph{BMJ}. 2017;356:i6460.
#' 
#' Lewis S, Clarke M. Forest plots: trying to see the wood and the trees. \emph{BMJ}. 2001; 322(7300):1479--80.
#' 
#' Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{BMJ}. 2011 342:d549--d549.
#' 
#' @seealso When a Bayesian meta-analysis was conducted, the prior and posterior distribution can be visualized using \code{\link{dplot.valmeta}}. 
#' Further, the running means and the presence of autocorrelation can be inspected using \code{\link{rmplot.valmeta}} and, respectively, 
#' \code{\link{acplot.valmeta}}.
#' 
#' @examples 
#' data(EuroSCORE)
#' fit <- valmeta(cstat=c.index, cstat.se=se.c.index, cstat.cilb=c.index.95CIl,
#'                cstat.ciub=c.index.95CIu, N=n, O=n.events, data=EuroSCORE)
#' plot(fit)
#' 
#' library(ggplot2)
#' plot(fit, theme=theme_grey())
#' 
#' @keywords meta-analysis discrimination calibration forest
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @method plot valmeta
#' @export
plot.valmeta <- function(x,  ...) {
  if (is.null(x$slab)) {
    yi.slab <- rownames(x$data)
  } else {
    yi.slab <- c(as.character(x$slab))
  }
  
  yi <- c(x$data[,"theta"])
  ci.lb <- c(x$data[,"theta.cilb"])
  ci.ub <- c(x$data[,"theta.ciub"])
  
  # Back-transform the raw data
  if (x$model == "normal/logit") {
    yi <- sapply(yi, inv.logit)
    ci.lb <- sapply(ci.lb, inv.logit)
    ci.ub <- sapply(ci.ub, inv.logit)
  } else if (x$model == "normal/log" | x$model == "poisson/log") {
    yi <- sapply(yi, exp)
    ci.lb <- sapply(ci.lb, exp)
    ci.ub <- sapply(ci.ub, exp)
  } 
  
  yi.ci <- cbind(ci.lb, ci.ub)
  
  forest(theta = yi, 
         theta.ci.lb = yi.ci[,"ci.lb"], 
         theta.ci.ub = yi.ci[,"ci.ub"], 
         theta.slab = yi.slab, 
         theta.summary = x$est, 
         theta.summary.ci.lb = x$ci.lb,
         theta.summary.ci.ub = x$ci.ub, 
         theta.summary.pi.lb = x$pi.lb,
         theta.summary.pi.ub = x$pi.ub,
         xlim = attr(x$data,'plot_lim'),
         refline = attr(x$data,'plot_refline'), 
         xlab = attr(x$data,'estimand'),  
         ...)
}



#' Plot the autocorrelation of a Bayesian meta-analysis model
#' 
#' Function to display autocorrelation of a fitted Bayesian meta-analysis model.
#' 
#' @param x An object of class \code{"valmeta"}
#' @param \ldots Additional arguments which are currently not used
#' @return A \code{ggplot} object.
#' 
#' @details 
#' Results are displayed for the estimated mean (\code{mu}) and standard-deviation (\code{tau}) of the meta-analysis model.
#' 
#' @examples 
#' \dontrun{
#' data(EuroSCORE)
#' 
#' fit <- valmeta(cstat=c.index, cstat.se=se.c.index, cstat.cilb=c.index.95CIl,
#'                cstat.ciub=c.index.95CIu, N=n, O=n.events, 
#'                data=EuroSCORE, method="BAYES", slab=Study)
#' acplot(fit)
#' } 
#' 
#' @keywords meta-analysis convergence autocorrelation
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @export
acplot.valmeta <- function(x, ...) {
  if (!("runjags" %in% class(x$fit))) {
    stop("The object 'x' does not represent a Bayesian analysis!")
  }
  if (!requireNamespace("ggmcmc", quietly = TRUE)) {
    stop("The package 'ggmcmc' is currently not installed!")
  } 
  
  P <- data.frame(
    Parameter = c("mu.tobs", "bsTau"),
    Label = c("mu", "tau"))
  
  acplot(x$fit$mcmc, P = P, greek = TRUE, ...)

}

#' Plot the running means of a Bayesian meta-analysis model
#' 
#' Function to display running means of a fitted Bayesian meta-analysis model.
#' 
#' @param x An object of class \code{"valmeta"}
#' @param \ldots Additional arguments which are currently not used
#' @return A \code{ggplot} object.
#' 
#' @details 
#' Results are displayed for the estimated mean (\code{mu}) and standard-deviation (\code{tau}) of the meta-analysis model.
#' 
#' @examples 
#' \dontrun{
#' data(EuroSCORE)
#' 
#' fit <- valmeta(cstat=c.index, cstat.se=se.c.index, cstat.cilb=c.index.95CIl,
#'                cstat.ciub=c.index.95CIu, N=n, O=n.events, 
#'                data=EuroSCORE, method="BAYES", slab=Study)
#' rmplot(fit)
#' } 
#' 
#' @keywords meta-analysis convergence
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @export
rmplot.valmeta <- function(x, ...) {
  if (!("runjags" %in% class(x$fit))) {
    stop("The object 'x' does not represent a Bayesian analysis!")
  }
  
  P <- data.frame(
    Parameter = c("mu.tobs", "bsTau"),
    Label = c("mu", "tau"))
  
  rmplot(x$fit$mcmc, P = P, greek = TRUE, ...)
}

#' Plot the prior and posterior distribution of a meta-analysis model
#' 
#' Function to generate plots for the prior and posterior distribution of a Bayesian meta-analysis.
#' 
#' @param x An object of class \code{"valmeta"}
#' @param par Character string to specify for which parameter a plot should be generated. Options are \code{"mu"} 
#' (mean of the random effects model) and \code{"tau"} (standard deviation of the random effects model).
#' @param distr_type Character string to specify whether the prior distribution (\code{"prior"}) or 
#' posterior distribution (\code{"posterior"}) should be displayed.
#' @param plot_type Character string to specify whether a density plot (\code{"dens"}) or 
#' histogram (\code{"hist"}) should be displayed.
#' @param \ldots Additional arguments which are currently not used
#' @return A \code{ggplot} object.
#' 
#' @examples 
#' \dontrun{
#' data(EuroSCORE)
#' 
#' # Meta-analysis of the concordance statistic
#' fit <- valmeta(cstat=c.index, cstat.se=se.c.index, cstat.cilb=c.index.95CIl,
#'                cstat.ciub=c.index.95CIu, N=n, O=n.events, 
#'                data=EuroSCORE, method="BAYES", slab=Study)
#' dplot(fit)
#' dplot(fit, distr_type = "posterior")
#' dplot(fit, par = "tau", distr_type = "prior")
#' 
#' # Meta-analysis of the O:E ratio
#' EuroSCORE.new <- EuroSCORE
#' EuroSCORE.new$n[c(1, 2, 5, 10, 20)] <-  NA
#' pars <- list(hp.tau.dist="dhalft",   # Prior for the between-study standard deviation
#'              hp.tau.sigma=1.5,       # Standard deviation for 'hp.tau.dist'
#'              hp.tau.df=3,            # Degrees of freedom for 'hp.tau.dist'
#'              hp.tau.max=10)          # Maximum value for the between-study standard deviation
#' fit2 <- valmeta(measure="OE", O=n.events, E=e.events, N=n, data=EuroSCORE.new,
#'                 method="BAYES", slab=Study, pars=pars)
#' dplot(fit2, plot_type = "hist")
#' } 
#' 
#' @keywords meta-analysis density distribution
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @export
dplot.valmeta <- function(x, par, distr_type, plot_type = "dens", ...) {
  if (!("runjags" %in% class(x$fit))) {
    stop("The object 'x' does not represent a Bayesian analysis!")
  }

  if (missing(par) & missing(distr_type)) {
    P <- data.frame(
      Parameter = c("prior_mu", "mu.tobs", "prior_bsTau", "bsTau"),
      Label = c("a) Prior distribution of the meta-analysis mean", 
                "b) Posterior distribution of the meta-analysis mean", 
                "c) Prior distribution of the between-study standard deviation",
                "d) Posterior distribution of the between-study standard deviation"))
  } else if (missing(par) & distr_type == "posterior") {
    P <- data.frame(
      Parameter = c("mu.tobs", "bsTau"),
      Label = c("Posterior distribution of the meta-analysis mean", "Posterior distribution of the between-study standard deviation"))
  } else if (missing(par) & distr_type == "prior") {
    P <- data.frame(
      Parameter = c("prior_mu", "prior_bsTau"),
      Label = c("Prior distribution of the meta-analysis mean", "Prior distribution of the between-study standard deviation"))
  } else if (par == "mu" & distr_type == "posterior") {
    P <- data.frame(
      Parameter = c("mu.tobs"),
      Label = c("Posterior distribution of the meta-analysis mean"))
  } else if (par == "mu" & distr_type == "prior") {
    P <- data.frame(
      Parameter = c("prior_mu"),
      Label = c("Prior distribution of the meta-analysis mean"))
  } else if (par == "tau") {
    P <- data.frame(
      Parameter = ifelse(distr_type == "posterior", "bsTau", "prior_bsTau"),
      Label = paste(ifelse(distr_type == "posterior", "Posterior", "Prior"), "distribution of the between-study standard deviation"))
  } else {
    stop("Invalid combination of 'par' and 'distr_type'")
  }
  
  dplot(x$fit$mcmc, P = P, plot_type = plot_type, ...)
}


#' Gelman-Rubin-Brooks plot
#' 
#' This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases. The code is adapted from
#' the R package coda.
#' 
#' @param x An mcmc object
#' @param confidence The coverage probability of the confidence interval for the potential scale reduction factor
#' @param \ldots Additional arguments which are currently not used
#' @return A \code{ggplot} object.
#' 
#' @examples 
#' \dontrun{
#' data(EuroSCORE)
#' 
#' # Meta-analysis of the concordance statistic
#' fit <- valmeta(cstat=c.index, cstat.se=se.c.index, cstat.cilb=c.index.95CIl,
#'                cstat.ciub=c.index.95CIu, N=n, O=n.events, 
#'                data=EuroSCORE, method="BAYES", slab=Study)
#' gelmanplot(fit)
#' } 
#' 
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @export
gelmanplot.valmeta <- function(x, confidence = 0.95, ...) {
  if (!("runjags" %in% class(x$fit))) {
    stop("The object 'x' does not represent a Bayesian analysis!")
  }
  
  P <- data.frame(
    Parameter = c("mu.tobs", "bsTau"),
    Label = c("mu", "tau"))
  
  return(gelmanplot(x$fit$mcmc, P = P, confidence = confidence, greek = TRUE, ...))
}