#Multivariate meta-analyse: http://blogs.sas.com/content/iml/2012/10/31/compute-the-log-determinant-of-a-matrix/ (directly take log-determinant)
#TODO: allow data transformations

#' Univariate meta-analysis
#' 
#' This function summarizes multiple estimates for a single parameter by assuming a fixed (i.e. common) 
#' effect or random effects across studies. The summary estimate is obtained by calculating a weighted mean that accounts for
#' sample size and (in case random effects are assumed) for between-study heterogeneity. 
#' 
#' @param r Vector of numerics containing the effect size of each study
#' @param r.se Vector of numerics containing the standard error of the effect sizes
#' @param r.vi Vector of numerics containing the sampling variance of the effect sizes
#' @param method Character string specifying whether a fixed-effect or a random-effects model should be fitted. 
#' A fixed-effect model is fitted when using \code{method="FE"}. Random-effects models are fitted by setting method equal 
#' to one of the following: \code{"REML"} (Default), \code{"DL"}, \code{"HE"}, \code{"SJ"}, \code{"ML"}, \code{"EB"}, 
#' \code{"HS"}, \code{"GENQ"} or \code{"BAYES"}. See 'Details'.
#' @param test Optional character string when \code{method!="BAYES"} to specify how test statistics and confidence intervals 
#' for the fixed effects should be computed. By default (\code{test="knha"}), the method by Knapp and Hartung (2003) is used for adjusting test statistics 
#' and confidence intervals.  Type '\code{?rma}' for more details.
#' @param labels Optional vector of characters containing the labels for the studies
#' @param na.action A function which indicates what should happen when the data contain NAs. 
#' Defaults to \code{"na.fail"}, other options are \code{"na.omit"}, \code{"na.exclude"} or \code{"na.pass"}.
#' @param n.chains Optional numeric specifying the number of chains to use in the Gibbs sampler (\code{method="BAYES"}). 
#' More chains will improve the sensitivity of the convergence diagnostic, but will cause the simulation to run more slowly. 
#' The default number of chains is 4.
#' @param pars Optional list with additional arguments. The width of confidence, credibility and prediction intervals is 
#' defined by \code{level} (defaults to 0.95). 
#' The following parameters configure the MCMC sampling procedure:  
#' \code{hp.mu.mean} (mean of the prior distribution of the random effects model, defaults to 0), 
#' \code{hp.mu.var} (variance of the prior distribution of the random effects model, defaults to 1E6), 
#' \code{hp.tau.min} (minimum value for the between-study standard deviation, defaults to 0), 
#' \code{hp.tau.max} (maximum value for the between-study standard deviation, defaults to 2), 
#' \code{hp.tau.sigma} (standard deviation of the prior distribution for the between-study standard-deviation), 
#' \code{hp.tau.dist} (prior distribution for the between-study standard-deviation. Defaults to \code{"dunif"}), 
#' \code{hp.tau.df} (degrees of freedom for the prior distribution for the between-study standard-deviation. 
#' Defaults to 3).
#' @param verbose If TRUE then messages generated during the fitting process will be displayed.
#' @param \dots Additional arguments that are passed to \pkg{rma} or \pkg{runjags} (if \code{method="BAYES"}).
#' 
#' @details Unless specified otherwise, all meta-analysis models assume random effects and are fitted  using restricted 
#' maximum likelihood estimation with the \pkg{metafor} package (Viechtbauer 2010).  Further, confidence intervals for 
#' the average performance are based on the Hartung-Knapp-Sidik-Jonkman method, to better account for the uncertainty 
#' in the estimated between-study heterogeneity (Debray 2016). A Bayesian meta-analysis can be performed by specifying 
#' \code{method="BAYES"}. In that case, the R packages \pkg{runjags} and \pkg{rjags} must be installed.]\cr
#' \cr
#' For random-effects models, a prediction interval for the pooled effect size is displayed. This interval predicts in what 
#' range future effect sizes will fall given what has already been observed (Higgins 2009, Riley 2011).  
#' 
#' \subsection{Bayesian meta-analysis models}{
#' For Bayesian meta-analysis models that involve the Gibbs sampler (\code{method="BAYES"}), the R packages \code{runjags} 
#' and \code{rjags} must be installed. The Bayesian approach uses an uninformative Normal prior for the mean and a 
#' uniform prior for the between-study variance of the pooled effect size (Higgins 2009). By default, the Normal prior 
#' has a mean of 0 and a variance of 1000. These hyperparameters can, however, be altered through the 
#' variables \code{hp.mu.mean} and \code{hp.mu.var} in the argument \code{pars}. The prior distribution of the between-study 
#' standard deviation is given by a uniform distribution, by default bounded between 0 and 100. 
#' }
#' 
#' @return An object of the class \code{uvmeta} for which many standard methods are available.
#' \describe{
##'  \item{"data"}{array with (transformed) data used for meta-analysis, and method(s) used for restoring missing information. }
##'  \item{"method"}{character string specifying the meta-analysis method.}
##'  \item{"est"}{estimated performance statistic of the model. For Bayesian meta-analysis, the posterior median is returned.}
##'  \item{"se"}{standard error (or posterior standard deviation) of the summary estimate.}
##'  \item{"tau2"}{estimated amount of (residual) heterogeneity. Always 0 when method=\code{"FE"}. For Bayesian meta-analysis, the posterior median is returned.}
##'  \item{"se.tau2"}{estimated standard error (or posterior standard deviation) of the between-study variation.}
##'  \item{"ci.lb"}{lower bound of the confidence (or credibility) interval of the summary estimate}
##'  \item{"ci.ub"}{upper bound of the confidence (or credibility) interval of the summary estimate}
##'  \item{"pi.lb"}{lower bound of the (approximate) prediction interval of the summary estimate}
##'  \item{"pi.ub"}{upper bound of the (approximate) prediction interval of the summary estimate}
##'  \item{"fit"}{the full results from the fitted model}
##'  \item{"slab"}{vector specifying the label of each study.}
##' }
#' 
#' @references
#' Biggerstaff BJ, Tweedie RL. Incorporating variability in estimates of heterogeneity in the random effects model 
#' in meta-analysis. \emph{Statistics in Medicine} 1997; \bold{16}: 753--768.
#' 
#' Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. A basic introduction to fixed-effect and random-effects 
#' models for meta-analysis. \emph{Research Synthesis Methods} 2010; \bold{1}: 97--111. \doi{10.1002/jrsm.12}
#' 
#' DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.
#' 
#' Graham PL, Moran JL. Robust meta-analytic conclusions mandate the provision of prediction intervals in 
#' meta-analysis summaries. \emph{Journal of Clinical Epidemiology} 2012; \bold{65}: 503--510.
#' 
#' Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine} 2002; 
#' \bold{21}: 1539--1558.
#' 
#' Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of random-effects meta-analysis. \emph{J R Stat Soc Ser A Stat Soc}. 
#' 2009;\bold{172}:137--59. \doi{10.1111/j.1467-985X.2008.00552.x}
#' 
#' Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. 
#' \emph{British Medical Journal} 2011; \bold{342}: d549. \doi{10.1136/bmj.d549}
#' 
#' Viechtbauer W. Conducting Meta-Analyses in R with the metafor Package. \emph{Journal of Statistical Software}. 
#' 2010; \bold{36}. \doi{10.18637/jss.v036.i03}
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @examples 
#' data(Roberts)
#' 
#' # Frequentist random-effects meta-analysis
#' fit1 <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts)))
#' summary(fit1)
#' plot(fit1) #show a forest plot
#' fit1
#' 
#' \dontrun{
#' # Bayesian random effects meta-analysis 
#' fit2 <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts), method="BAYES"))
#' plot(fit2)
#' }
#' 
#' @keywords univariate fixed-effect random-effects meta-analysis heterogeneity
#' 
#' @import metafor
#' @export

uvmeta <- function(r, r.se, r.vi, method="REML", test="knha", labels, na.action, 
                   n.chains=4, pars, verbose=FALSE, ...) 
  UseMethod("uvmeta")

#' @export
uvmeta.default <- function(r, r.se, r.vi, method="REML", test="knha", labels, na.action, 
                           n.chains=4, pars, verbose=FALSE, ...)
{
  out <- list()
  out$call <- match.call()
  out$method <- method
  out$test <- test
  class(out) <- "uvmeta"
  
  pars.default <- .initiateDefaultPars(pars)

  
  # Check if we need to load runjags
  if (method=="BAYES") {
    if (!requireNamespace("runjags", quietly = TRUE)) {
      stop("The package 'runjags' is currently not installed!")
    } 
    if (!requireNamespace("rjags", quietly = TRUE)) {
      stop("The package 'rjags' is currently not installed!")
    } 
    if (n.chains<1 | n.chains%%1!=0) {
      stop("Invalid number of chains specified for the Gibbs sampler!")
    }
    if (pars.default$hp.tau.min < 0) {
      stop("Invalid value for hyperparameter 'hp.tau.min")
    }
    if (pars.default$hp.tau.max < pars.default$hp.tau.min) {
      stop("Invalid value for hyperparameter 'hp.tau.max")
    }
  }
  
  if (!missing(r.se)) {
    if (length(r)!=length(r.se)) {
      stop("The vectors 'r' and 'r.se' have different lengths!")
    } 
  } else if (!missing(r.vi)) {
    if (length(r)!=length(r.vi)) {
      stop("The vectors 'r' and 'r.vi' have different lengths!")
    }
    r.se <- sqrt(r.vi)
  }
  
  ds <- as.data.frame(cbind(as.vector(r),as.vector(r.se)))
  colnames(ds) <- c("theta","theta.se")
  
  if (!missing(labels)) {
    if (length(labels) != length(r))
      stop("The vectors 'labels' and 'r' have different lengths!")
    
    out$slab <- make.unique(as.character(labels))
  } else {
    out$slab <- paste("Study",seq(1, length(r)))
  }
  rownames(ds) = out$slab
  
  if (missing(na.action)) 
    na.action <- "na.fail"
  if (length(na.action)) 
    ds <- do.call(na.action, list(ds))
  
  
  
  # Define quantiles for calculating intervals
  out$level <- pars.default$level
  quantiles <- c((1-pars.default$level)/2, 0.50, (1-((1-pars.default$level)/2)))
  
  if (out$level < 0 | out$level > 1) {
    stop("Invalid value for 'level'!")
  } 
  
  
  #############################################################################
  # Start analyses
  #############################################################################
  numstudies <- dim(ds)[1]

  if (numstudies < 3) {
    warning("There are very few primary studies!")
  }
  
  if (method == "BAYES") {
      return(run_Bayesian_REMA(call = match.call(),
                               method = method,
                               data = ds, 
                               pars = pars.default, 
                               FUN_generate_bugs = .generateBugsREMA,
                               n.chains = n.chains, 
                               verbose = verbose, ...))
  }
  
  
  if (method != "BAYES") { 
    fit <- metafor::rma(yi=r, sei=r.se, method=method, test=test, slab=out$slab, ...) 
    preds <- predict(fit, level=pars.default$level)
    
    # We don't use the prediction intervals from metafor, as they are based on a Normal distribution
    predint <- calcPredInt(coefficients(fit), sigma2=fit$se**2, tau2=fit$tau2, k=fit$k, level=pars.default$level)
    
    out$est <- as.numeric(coefficients(fit))
    out$se  <- fit$se
    out$tau2 <- fit$tau2
    out$se.tau2 <- fit$se.tau2
    out$ci.lb <- preds$ci.lb
    out$ci.ub <- preds$ci.ub
    out$pi.lb <- ifelse(method=="FE", preds$ci.lb, predint$lower)
    out$pi.ub <- ifelse(method=="FE", preds$ci.ub, predint$upper)

    out$fit <- fit
    out$numstudies <- fit$k
    
  }
  #attr(out$results,"level") <- pars.default$level
  out$data <- ds
  out$numstudies <- dim(ds)[1]
  out$na.action <- na.action
  class(out) <- "uvmeta"
  
  return(out)
}

#' Forest Plots
#' 
#' Function to create forest plots for objects of class \code{"uvmeta"}.
#' 
#' @param x An object of class \code{"uvmeta"}
#' @param sort By default, studies are ordered by ascending effect size (\code{sort="asc"}). For study ordering by descending
#' effect size, choose \code{sort="desc"}. For any other value, study ordering is ignored.
#' @param \dots Additional arguments which are passed to \link{forest}.
#' 
#' @details The forest plot shows the performance estimates of each validation with corresponding confidence 
#' intervals. A polygon is added to the bottom of the forest plot, showing the summary estimate based on the model. 
#' A 95\% prediction interval is added by default for random-effects models,  the dotted line indicates its (approximate) bounds
#' 
#' @note Full lines indicate confidence intervals or credibility intervals (in case of a Bayesian meta-analysis). Dashed
#' lines indicate prediction intervals. The width of all intervals is defined by the significance level chosen during 
#' meta-analysis. 
#' 
#' @references \itemize{
#' \item Lewis S, Clarke M. Forest plots: trying to see the wood and the trees. \emph{BMJ}. 2001; 322(7300):1479--80.
#' \item Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{BMJ}. 2011 342:d549--d549.
#' }
#' 
#' @examples 
#' data(Roberts)
#' 
#' # Frequentist random-effects meta-analysis
#' fit <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts)))
#' plot(fit) 
#' 
#' @keywords meta-analysis forest
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @importFrom graphics par plot axis polygon points lines abline
#' @import ggplot2
#' 
#' @method plot uvmeta
#' @export
plot.uvmeta <- function(x, sort="asc", ...) {
  level <- x$level
  quantiles <- c((1-level)/2, (1-((1-level)/2)))
  
  # Reconstruct confidence intervals for study effects with the level used during calculation of summary effect
  yi.ci <- x$data[,"theta"]+t(qnorm(quantiles)*matrix(rep(x$data[,"theta.se"],length(quantiles)),nrow=(length(quantiles)), ncol=dim(x$data)[1],byrow=T))
  
  #Extract data
  yi <- c(x$data[,"theta"])

  if (is.null(x$slab)) {
    yi.slab <- rownames(x$data)
  } else {
    yi.slab <- c(as.character(x$slab))
  }
  
  forest(theta = yi, 
         theta.ci.lb = yi.ci[,1], 
         theta.ci.ub = yi.ci[,2],
         theta.slab = yi.slab, 
         theta.summary = x$est, 
         theta.summary.ci.lb = x$ci.lb,
         theta.summary.ci.ub = x$ci.ub,
         theta.summary.pi.lb = x$pi.lb,
         theta.summary.pi.ub = x$pi.ub,
         sort = sort, ...)
}


#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method print uvmeta
#' @export
print.uvmeta <- function(x, ...)
{
  text.model <- if (x$method=="FE") "Fixed" else "Random"
  text.ci <- if(x$method=="BAYES") "credibility" else "confidence"
  text.pi <- if(x$method=="BAYES") "" else "(approximate)"
  
  if (x$method!="FE") {
    cat(paste("Summary estimate with ", x$level*100, "% ", text.ci, " and ", text.pi, " ", 
              x$level*100, "% prediction interval:\n\n", sep=""))
    results=c(Estimate=x$est, CIl=x$ci.lb, CIu=x$ci.ub, PIl=x$pi.lb, PIu=x$pi.ub)
  } else {
    cat(paste("Summary estimate with ", x$level*100, "% ", text.ci, ":\n\n", sep=""))
    results=c(Estimate=x$est, CIl=x$ci.lb, CIu=x$ci.ub)
  }
  print(results)
}

#' Summarizing Univariate Meta-Analysis Models
#' 
#' This function provides summary estimates of a fitted univariate meta-analysis model.
#' 
#' @method summary uvmeta
#' 
#' @param  object  An object of class \code{"uvmeta"}
#' @param \dots Optional arguments to be passed on to other functions
#' 
#' @references  \itemize{
#' \item Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. A basic introduction to fixed-effect and random-effects models for meta-analysis. \emph{Research Synthesis Methods} 2010; \bold{1}: 97--111.
#' \item DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.
#' \item Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{British Medical Journal} 2011; \bold{342}: d549.
#' }
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' @seealso \code{\link{uvmeta}}
#' 
#' @export
#' @keywords DerSimonian  Laird univariate random-effects meta-analysis 
summary.uvmeta <- function(object, ...)
{
    cat("Call:\n")
    print(object$call)
    if (object$method=="FE")  cat(paste("\nFixed effects summary:\t",round(object$est,5))," (SE: ",round(object$se,5), ")",sep="")
    if (object$method!="FE") {
        cat(paste("\nRandom effects summary:\t",round(object$est,5))," (SE: ",round(object$se,5), ")",sep="")
        cat(paste("\nTau squared: \t\t",round(object$tau2,5))," (SE: ",round(object$se.tau2,5), ")",sep="")
    }
}

#' Plot the prior and posterior distribution of a meta-analysis model
#' 
#' Function to generate plots for the prior and posterior distribution of a Bayesian meta-analysis.
#' 
#' @param x An object of class \code{"uvmeta"}
#' @param par Character string to specify for which parameter a plot should be generated. Options are \code{"mu"} 
#' (mean of the random effects model) and \code{"tau"} (standard deviation of the random effects model).
#' @param distr_type Character string to specify whether the prior distribution (\code{"prior"}) or 
#' posterior distribution (\code{"posterior"}) should be displayed.
#' @param plot_type Character string to specify whether a density plot (\code{"dens"}) or 
#' histogram (\code{"hist"}) should be displayed.
#' @param \ldots Additional arguments which are currently not used
#' 
#' @examples 
#' \dontrun{
#' data(Roberts)
#' 
#' fit <- with(Roberts, uvmeta(r=SDM, r.se=SE, method="BAYES"))
#' 
#' dplot(fit)
#' dplot(fit, distr_type = "posterior")
#' dplot(fit, par = "tau", distr_type = "prior")
#' dplot(fit, plot_type = "hist")
#' } 
#' 
#' @keywords meta-analysis density distribution
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @export
dplot.uvmeta <- function(x, par, distr_type, plot_type = "dens", ...) {
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
#' data(Roberts)
#' 
#' fit <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts), method="BAYES"))
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
acplot.uvmeta <- function(x, ...) {
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
#' data(Roberts)
#' 
#' fit <- with(Roberts, uvmeta(r=SDM, r.se=SE, labels=rownames(Roberts), method="BAYES"))
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
rmplot.uvmeta <- function(x, ...) {
  if (!("runjags" %in% class(x$fit))) {
    stop("The object 'x' does not represent a Bayesian analysis!")
  }
  
  P <- data.frame(
    Parameter = c("mu.tobs", "bsTau"),
    Label = c("mu", "tau"))
  
  rmplot(x$fit$mcmc, P = P, greek = TRUE, ...)
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
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @export
gelmanplot.uvmeta <- function(x, confidence = 0.95, ...) {
  if (!("runjags" %in% class(x$fit))) {
    stop("The object 'x' does not represent a Bayesian analysis!")
  }
  
  P <- data.frame(
    Parameter = c("mu.tobs", "bsTau"),
    Label = c("mu", "tau"))
  
  return(gelmanplot(x$fit$mcmc, P = P, confidence = confidence, greek = TRUE, ...))
}