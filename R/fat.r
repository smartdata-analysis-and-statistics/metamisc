### To add / change:
# Add plot support for P-FPV, D-FIV and D-FAV
# Add subsections for Details
# 
#' Regression tests for detecting funnel plot asymmetry
#'
#' The presence of small-study effects is a common threat to systematic reviews and meta-analyses, especially when 
#' it is due to publication bias, which occurs when small primary studies are more likely to be reported (published) 
#' if their findings were positive. The presence of small-study effects can be verified by visual inspection of 
#' the funnel plot, where for each included study of the meta-analysis, the estimate of the reported effect size is 
#' depicted against a measure of precision or sample size. 
#' The premise is that the scatter of plots should reflect a funnel shape, if small-study 
#' effects do not exist. However, when small studies are predominately in one direction (usually the 
#' direction of larger effect sizes), asymmetry will ensue.\cr \cr
#' The \code{\link{fat}} function implements several tests for detecting funnel plot asymmetry, 
#' which can be used when the presence of between-study heterogeneity in treatment effect is relatively low.
#' 
#' @param b Vector with the effect size of each study. Examples are log odds ratio, log hazards ratio, 
#' log relative risk. 
#' @param b.se Optional vector with the standard error of the effect size of each study
#' @param n.total Optional vector with the total sample size of each study
#' @param d.total Optional vector with the total number of observed events for each study
#' @param d1 Optional vector with the total number of observed events in the exposed groups
#' @param d2 Optional vector with the total number of observed events in the unexposed groups
#' @param method Method for testing funnel plot asymmetry, defaults to \code{"E-FIV"} (Egger's test with 
#' multiplicative dispersion). Other options are \code{E-UW}, \code{M-FIV}, \code{M-FPV}, \code{D-FIV} and 
#' \code{D-FAV}. More info in "Details"
#'
#' @details 
#' 
#' \subsection{Egger regression method}{A common method to test the presence of small-study effects is given as the 
#' following unweighted regression model (\code{method="E-UW"}, Egger 1997): 
#' \deqn{\hat{b}_k = \beta_0 + \beta_1\, \widehat \mathrm{SE}(\hat{b}_k) + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}\left(0, \sigma^2 \right) }{b = B0 + B1*b.se + e;  e~N(0, s^2)}
#' Whereas \eqn{\beta_0}{B0} indicates the size and direction of the treatment effect, \eqn{\beta_1}{B1} provides 
#' a measure of asymmetry; the larger its deviation from zero the more pronounced the asymmetry. Otherwise, if 
#' \eqn{\beta_1=0}{B1=0}, there is no association between the estimated effect sizes \eqn{\hat{b}_k}{b} and their 
#' corresponding estimates for the standard error \eqn{\widehat \mathrm{SE}(\hat{b}_k)}{b.se} among the reported 
#' studies, indicating no asymmetry and thus no small-study effects. \cr \cr
#' It is possible to allow for between-study heterogeneity by adopting a multiplicative overdispersion parameter 
#' by which the variance in each study is multiplied (\code{method="E-FIV"}, Sterne 2000):
#' \deqn{\hat{\beta}_k = a + b\, \widehat \mathrm{SE}(\hat{\beta}_k) + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}(0, \phi \; \widehat \mathrm{var}(\hat{\beta}_k))}{b = B0 + B1*b.se + e;  e~N(0, P*b.se^2)}
#' Unfortunately, both tests are known to be intrinsically biased because: (i) the independent variable is subject to sampling 
#' variability; (ii) the standardized treatment effect is correlated with its estimated precision; and 
#' (iii) for binary data, the independent regression variable is a biased estimate of the true precision, 
#' with larger bias for smaller sample sizes (Macaskill 2001).
#' }
#' 
#' \subsection{Macaskill regression method}{
#' To overcome the problems with the Egger approach, Macaskill et al. consider fitting a regression directly
#' to the data using the treatment effect as the dependent variable, and study size (\eqn{n_k}{n.total}) as the 
#' independent variable. Again, the observations are weighted by the inverse variance of the estimate
#' to allow for possible heteroscedasticity  (\code{method="M-FIV"}, Macaskill 2001):
#' \deqn{\hat{\beta}_k = a + b \,n_k + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}(0, \phi \; \widehat \mathrm{var}(\hat{\beta}_k))}{b = B0 + B1*n.total + e;  e~N(0, P*b.se^2)}
#' Macaskill et al. also proposed an alternative test where a 'pooled' estimate of the outcome proportion is used
#' for the variance \eqn{\widehat \mathrm{var}(\hat{b}_k)}{b.se^2} (\code{method="M-FPV"}, Macaskill 2001):
#' \deqn{\hat{\beta}_k = a + b \,n_k + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}\left(0, \phi \; \frac{1}{d_k (1-d_k/n_k)}\right)}{b = B0 + B1*n.total + e;  e~N(0, P/(d.total * (1-d.total/n.total)))}
#' For studies with zero events, a continuity correction is applied by adding 0.5 to all cell counts.
#' }
#' 
#' \subsection{Peters regression method}{
#' A modification of Macaskill's test was proposed by Peters et al. to obtain more balanced type-I error rates 
#' in the tail probability areas  (\code{method="P-FPV"}, Peters 2006):
#' \deqn{\hat{\beta}_k = a + b \,\frac{1}{n_k} + \epsilon_k \;,\; \epsilon_k \sim \mathcal{N}\left(0, \phi \; \frac{1}{d_k (1-d_k/n_k)}\right)}{b = B0 + B1/n.total + e;  e~N(0, P/(d.total * (1-d.total/n.total)))}
#' Again, 0.5 is added to all cells for studies with zero events.
#' }
#' 
#' \subsection{Debray regression method}{
#' Because the use of aforementioned tests may be less appropriate in the presence of survival data, Debray et al. 
#' proposed using the total number of events (\eqn{d_k}{d.total}) as independent variable (\code{method="D-FIV"}, Debray 2017):
#' \deqn{\hat{\beta}_k = a + b\, \frac{1}{d_k} + \epsilon_k  \;,\; \epsilon_k \sim \mathcal{N}(0, \phi \; \widehat \mathrm{var}(\hat{\beta}_k))}{b = B0 + B1/d.total + e;  e~N(0, P*b.se^2)}
#' For studies with zero events, the total number of observed events is set to 1.
#' Alternatively, when \eqn{\widehat \mathrm{var}(\hat{\beta}_k)}{b.se} is unknown or derived from small samples, 
#' Debray at al.proposed to use the following regression model (\code{method="D-FAV"}, Debray 2017):
#' \deqn{\hat{\beta}_k = a + b\, \frac{1}{d_k} + \epsilon_k  \;,\; \epsilon_k \sim \mathcal{N}\left(0, \phi \; \left(\frac{1}{d_{k1}}+\frac{1}{d_{k2}}\right)\right)}{b = B0 + B1/d.total + e;  e~N(0, P/(1/d1 + 1/d2))}
#' }
#' 
#' @return a list containing the following entries:
#' \describe{
##'  \item{"pval"}{A two-sided P-value indicating statistical significance of the funnel plot asymettry test. 
##'  Values below the significance level (usually defined as 10\%) support the presence of funnel plot asymmetry,
##'  and thus small-study effects.  }
##'  \item{"model"}{A fitted \code{glm} object, representing the estimated regression model used for testing funnel
##'  plot asymmetry.}
##' }
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @references Debray TPA, Moons KGM, Riley RD. Detecting small-study effects and funnel plot asymmetry in meta-analysis of 
#' survival data: a comparison of new and existing tests. Res Syn Meth. 2018;9(1):41--50.\cr
#' \cr
#' Egger M, Davey Smith G, Schneider M, Minder C. Bias in meta-analysis detected by a simple, graphical test. 
#' \emph{BMJ}. 1997;315(7109):629--34. \cr
#' \cr
#' Macaskill P, Walter SD, Irwig L. A comparison of methods to detect publication bias in meta-analysis. 
#' \emph{Stat Med}. 2001;20(4):641--54.\cr 
#' \cr
#' Peters JL, Sutton AJ, Jones DR, Abrams KR, Rushton L. Comparison of two methods to detect publication bias 
#' in meta-analysis. \emph{JAMA}. 2006 Feb 8;295(6):676--80.\cr
#' \cr 
#' Sterne JA, Gavaghan D, Egger M. Publication and related bias in meta-analysis: power of statistical tests 
#' and prevalence in the literature. \emph{J Clin Epidemiol}. 2000;53(11):1119--29. 
#' 
#' @seealso \code{\link{plot.fat}}
#'
#' @examples 
#' data(Fibrinogen)
#' b <- log(Fibrinogen$HR)
#' b.se <- ((log(Fibrinogen$HR.975) - log(Fibrinogen$HR.025))/(2*qnorm(0.975)))
#' n.total <- Fibrinogen$N.total
#' d.total <- Fibrinogen$N.events
#' 
#' fat(b=b, b.se=b.se)
#' fat(b=b, b.se=b.se, d.total=d.total, method="D-FIV")
#' 
#' # Note that many tests are also available via metafor
#' require(metafor)
#' fat(b=b, b.se=b.se, n.total=n.total, method="M-FIV")
#' regtest(x=b, sei=b.se, ni=n.total, model="lm", predictor="ni") 
#'
#' @import stats
#' @importFrom stats pt qnorm
#' @importFrom metafor rma
#' @importFrom plyr round_any
#' 
#' @export
fat <- function(b, b.se, n.total, d.total, d1, d2, method="E-FIV") 
{
  if (missing(b)) {
    stop ("No values given for 'b'")
  }
  
  # Identify studies with complete information
  if (method == "E-UW") {
    if (missing(b.se)) {
      stop ("No values given for 'b.se'")
    }
    if (length(b) != length(b.se)) {
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(b.se))
    ds <- data.frame("y" = b, 
                     "x" = b.se
                     )
  } else if (method== "E-FIV") {
    if (missing(b.se)) {
      stop ("No values given for 'b.se'")
    }
    if (length(b) != length(b.se)) {
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(b.se))
    ds <- data.frame("y" = b, 
                     "x" = b.se, 
                     "w" = (1/(b.se**2))
                     )
  } else if (method == "M-FIV") {
    if (missing(b.se)) {
      stop ("No values given for 'b.se'")
    }
    if (missing(n.total)) {
      stop ("No values given for 'n.total'")
    }
    if (length(b) != length(b.se)) {
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    }
    if (length(b) != length(n.total)) {
      stop("Incompatible vector sizes for 'b' and 'n.total'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(b.se) & !is.na(n.total))
    ds <- data.frame("y" = b, 
                     "x" = n.total, 
                     "w" = (1/(b.se**2))
                     )
  } else if (method=="M-FPV") {
    if (missing(n.total)) {
      stop ("No values given for 'n.total'")
    }
    if (missing(d.total)) {
      stop ("No values given for 'd.total'")
    }
    if (length(b) != length(n.total)) {
      stop("Incompatible vector sizes for 'b' and 'n.total'!")
    }
    if (length(b) != length(d.total)) {
      stop("Incompatible vector sizes for 'b' and 'd.total'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(d.total) & !is.na(n.total))
    
    # Consider continuity corrections
    d.total.cc <- d.total
    d.total.cc[d.total==0] <- 1 #0.5 event in exposed group and 0.5 event in non-exposed group
    n.total[d.total==0] <- n.total[d.total==0]+2 #2*0.5 in the events, and 2*0.5 in the non-events
    
    ds <- as.data.frame(cbind(b, n.total, (d.total.cc*(1-d.total.cc/n.total))))
    colnames(ds) <- c("y","x","w")
  } else if (method=="P-FPV") {
    if (missing(n.total)) {
      stop ("No values given for 'n.total'")
    }
    if (missing(d.total)) {
      stop ("No values given for 'd.total'")
    }
    if (length(b) != length(n.total)) {
      stop("Incompatible vector sizes for 'b' and 'n.total'!")
    }
    if (length(b) != length(d.total)) {
      stop("Incompatible vector sizes for 'b' and 'd.total'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(d.total) & !is.na(n.total))
    
    # Consider continuity corrections
    d.total.cc <- d.total
    d.total.cc[d.total==0] <- 1 #0.5 event in exposed group and 0.5 event in non-exposed group
    n.total[d.total==0] <- n.total[d.total==0]+2 #2*0.5 in the events, and 2*0.5 in the non-events
    
    ds <- data.frame("y" = b, 
                     "x" = 1/n.total, 
                     "w" = (d.total.cc*(1-d.total.cc/n.total))
                     )
  } else if (method=="D-FIV") {
    if (missing(b.se)) {
      stop ("No values given for 'b.se'")
    }
    if (missing(d.total)) {
      stop ("No values given for 'd.total'")
    }
    if (length(b) != length(b.se)) {
      stop("Incompatible vector sizes for 'b' and 'b.se'!")
    }
    if (length(b) != length(d.total)) {
      stop("Incompatible vector sizes for 'b' and 'd.total'!")
    }
    studies.complete <- c(!is.na(b) & !is.na(b.se) & !is.na(d.total))
    
    # Consider continuity corrections
    d.total.cc <- d.total
    d.total.cc[d.total==0] <- 1 #0.5 event in exposed group and 0.5 event in non-exposed group

    ds <- data.frame("y" = b, 
                     "x" = 1/d.total.cc, 
                     "w" = (1/(b.se**2))
                     )
  } else if (method=="D-FAV") {
    if (missing(d1)) {
      stop ("No values given for 'd1'")
    }
    if (missing(d2)) {
      stop ("No values given for 'd2'")
    }
    if (length(b) != length(d1)) {
      stop("Incompatible vector sizes for 'b' and 'd1'!")
    }
    if (length(b) != length(d2)) {
      stop("Incompatible vector sizes for 'b' and 'd2'!")
    }
    if (!missing(d.total)) {
      if (sum(d1+d2!=d.total) > 0)
        stop("Incompatible information between 'd.total', 'd1' and 'd2'")
    }
    studies.complete <- c(!is.na(b) & !is.na(d1) & !is.na(d2))
    
    # Consider continuity corrections
    d1.cc <- d1
    d2.cc <- d2
    d1.cc[(d1==0 | d2==0)] <- d1.cc[(d1==0 | d2==0)]+0.5 #0.5 event in exposed group and 0.5 event in non-exposed group
    d2.cc[(d1==0 | d2==0)] <- d2.cc[(d1==0 | d2==0)]+0.5
    
    ds <- data.frame("y" = b, 
                     "x" = 1/(d1.cc+d2.cc),  
                     "w" = 1/((1/d1.cc)+(1/d2.cc))
                     )
  } 
  else {
    stop("Method for testing funnel plot asymmetry not supported")
  }
  
  # Identify which studies can be used
  nstudies <- sum(studies.complete)
  
  # Omit sudies with missing information  
  ds <- ds[studies.complete,]
  
  if (nstudies < length(studies.complete)) {
    warning("Some studies were removed due to missing data!")
  }
  
  # Get the fixed effect summary estimate
  res <- NULL
  if (!missing(b.se)) {
    res <- rma(yi = b[studies.complete], sei = b.se[studies.complete], method = "FE")
  }
  

  
  if (method %in% c("E-FIV", "M-FIV", "M-FPV", "P-FPV", "D-FIV", "D-FAV")) {
    suppressWarnings(m.fat <- try(glm(y~x, weights=ds$w, data=ds), silent=T))
  } else if (method=="E-UW")  {
    suppressWarnings(m.fat <- try(glm(y~x, data=ds), silent=T))
  } else {
    stop("Method for testing funnel plot asymmetry currently not implemented")
  }
  
  if ("try-error" %in% attr(m.fat,"class")) {
    warning("Estimation of the regression model unsuccessful, P-value omitted.")
    t.fat <- NA
    p.fat <- NA
  } else {
    t.fat <- coefficients(m.fat)[2]/sqrt(diag(vcov(m.fat))[2])
    p.fat <- 2*pt(-abs(t.fat),df=(nstudies-2))
  }

  out <- list()
  out$call <- match.call()
  out$method <- method
  out$tval <- t.fat
  out$pval <- p.fat
  out$fema <- res
  out$df <- nstudies-2
  out$model <- m.fat
  class(out) <- "fat"
  return(out)
}

#' @method print fat
#' @export
print.fat <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Call: ");                       
  print(x$call); 
  if (!is.null(x$fema)) {
    cat(c("\nFixed effect summary estimate: ", round(x$fema$b, digits = digits), " \n"))
  }
  cat("\n")
  cat(paste("test for funnel plot asymmetry: t =", round(x$tval, digits = digits), ", df = ", x$df, ", ", sep=""))
  cat(paste("p = ", round(x$pval, digits = digits), "\n", sep=""))
  
}

#' Display results from the funnel plot asymmetry test
#' 
#' Generates a funnel plot for a fitted \code{fat} object.
#' @param x An object of class \code{fat}
#' @param ref A numeric value indicating the fixed or random effects summary estimate. If no value is provided
#' then it will be retrieved from a fixed effects meta-analysis (if possible).
#' @param confint A logical indicator. If \code{TRUE}, a confidence interval will be displayed for the estimated
#' regression model (based on a Student-T distribution)
#' @param confint.level Significance level for constructing the confidence interval.
#' @param confint.alpha A numeric value between 0 and 1 indicating the opacity for the confidence region.
#' @param confint.col The color for filling the confidence interval. Choose \code{NA} to leave polygons unfilled. 
#' If \code{confint.density} is specified with a positive value this gives the color of the shading lines. 
#' @param confint.density The density of shading lines, in lines per inch. The default value of \code{NULL} means 
#' that no shading lines are drawn. A zero value of density means no shading nor filling whereas negative values 
#' and \code{NA} suppress shading (and so allow color filling).
#' @param xlab A title for the x axis
#' @param add.pval Logical to indicate whether a P-value should be added to the plot
#' @param ... Additional arguments. 
#' 
#' @examples 
#' data(Fibrinogen)
#' b <- log(Fibrinogen$HR)
#' b.se <- ((log(Fibrinogen$HR.975) - log(Fibrinogen$HR.025))/(2*qnorm(0.975)))
#' n.total <- Fibrinogen$N.total
#' 
#' # A very simple funnel plot
#' plot(fat(b=b, b.se=b.se), xlab = "Log hazard ratio")
#' 
#' # Plot the funnel for an alternative test
#' plot(fat(b=b, b.se=b.se, n.total=n.total, method="M-FIV"), xlab = "Log hazard ratio")
#' 
#' 
#' @import ggplot2
#' @importFrom stats qt
#' @importFrom graphics plot axis polygon points lines box abline
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' @author Frantisek Bartos <f.bartos96@gmail.com>
#' @method plot fat
#' 
#' @export
plot.fat <- function(x, ref, confint = TRUE, confint.level = 0.10, confint.col = "skyblue", confint.alpha = .50,
                     confint.density = NULL,
                     xlab = "Effect size", add.pval = TRUE, ...) {
  if (!inherits(x, "fat"))
    stop("Argument 'x' must be an object of class \"fat\".")
  if (confint.level < 0 | confint.level > 1) {
    stop("Argument 'confint.level' must be between 0 and 1.")
  }
  if (confint.alpha < 0 | confint.alpha > 1) {
    stop("Argument 'confint.alpha' must be between 0 and 1.")
  }
  
  y <- NULL
  
  xval <- x$model$data[, "y"]
  if (x$method %in% c("E-UW", "E-FIV")) {
    ylab <- "Standard error"
    yval <- (x$model$data[, "x"])
    ylim <- rev(c(0, max(yval, na.rm = T)))
    xlim <- c(min(c(0, xval)), max(xval))
  } else if (x$method %in% c("M-FIV")) {
    ylab <- "Sample size"
    yval <- (x$model$data[, "x"])
    ylim <- (c(0, max(yval, na.rm = T)))
    xlim <- c(min(c(0, xval)), max(xval))
  } else if (x$method == "P-FPV") {
    ylab <- "Sample size"
    yval <- (x$model$data[, "x"])
    ylim <- rev(c(0, max(yval, na.rm = T)))
    xlim <- c(min(c(0, xval)), max(xval))
    step <- ((max(yval) - min(yval))/5)
    yax  <- c(plyr::round_any(1/min(yval), 10^(sapply(round(1/min(yval)), nchar) - 1)),
              plyr::round_any(1/seq(step, 4 * step, by = step), 10), plyr::round_any(1/max(yval), 10))
  } else if (x$method == "D-FIV") {
    ylab <- "Total events"
    yval <- (x$model$data[, "x"])
    ylim <- rev(c(0, max(yval, na.rm = T)))
    xlim <- c(min(c(0, xval)), max(xval))
    step <- ((max(yval) - min(yval))/4)
    yax  <- c(plyr::round_any(1/min(yval), 10^(sapply(round(1/min(yval)), nchar) - 1)),
              plyr::round_any(1/seq(step, 4 * step, by = step), 10), plyr::round_any(1/max(yval), 10))
  } else if (x$method == "D-FAV") {
    ylab <- "Total events"
    yval <- (x$model$data[, "x"])
    ylim <- rev(c(0, max(yval, na.rm = T)))
    xlim <- c(min(c(0, xval)), max(xval))
    step <- ((max(yval) - min(yval))/4)
    yax  <- c(plyr::round_any(1/min(yval),10^(sapply(round(1/min(yval)), nchar) - 1)),
              plyr::round_any(1/seq(step, 4 * step, by = step), 10), plyr::round_any(1/max(yval), 10))
  } else {
    stop("Plot not supported!")
  }
  
  newdata <- sort(c(-max(x$model$data[, "x"]), x$model$data[,"x"], 2 * max(x$model$data[, "x"])))
  newdata <- as.data.frame(cbind(seq(min(newdata), max(newdata), length.out = 500), NA))
  colnames(newdata) <- c("x", "y")
  predy <- predict(x$model, newdata = newdata, se.fit = T)
  predy.mean <- predy$fit
  predy.lowerInt <- as.vector(predy$fit + qt(confint.level/2,  df = x$df) * predy$se.fit)
  predy.upperInt <- as.vector(predy$fit + qt((1 - confint.level/2),  df = x$df) * predy$se.fit)
  
  # restricting plotting range to the selection
  predy.upperInt[predy.upperInt < min(pretty(range(xlim)))] <- min(pretty(range(xlim)))
  predy.upperInt[predy.upperInt > max(pretty(range(xlim)))] <- max(pretty(range(xlim)))
  predy.lowerInt[predy.lowerInt < min(pretty(range(xlim)))] <- min(pretty(range(xlim)))
  predy.lowerInt[predy.lowerInt > max(pretty(range(xlim)))] <- max(pretty(range(xlim)))
  newdata[, "x"][newdata[, "x"] < min(pretty(range(ylim)))] <- min(pretty(range(ylim)))
  newdata[, "x"][newdata[, "x"] > max(pretty(range(ylim)))] <- max(pretty(range(ylim)))
  
  p <- ggplot2::ggplot(data = data.frame(x = xval, y = yval))
  
  if (confint) {
    
    
    p <- p + ggplot2::geom_polygon(
      mapping = ggplot2::aes(
        x = x,
        y = y
      ),
      data = data.frame(
        x = c(
          predy.upperInt,
          rev(predy.lowerInt)),
        y = c(
          newdata[, "x"],
          rev(newdata[, "x"]))
      ),
      fill  = confint.col,
      alpha = confint.alpha
    )
  }
  
  p <- p +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = x, y = y),
      shape   = 19
    ) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(x = x, y = y),
      data    = data.frame(
        x = predy.mean[newdata[, "x"] > min(pretty(range(ylim))) & newdata[, "x"] < max(pretty(range(ylim)))],     # another plotting range restriction
        y = newdata[, "x"][newdata[, "x"] > min(pretty(range(ylim))) & newdata[, "x"] < max(pretty(range(ylim)))]  # another plotting range restriction  
      ),
      linetype = 2)
  
  
  
  if (missing(ref)) {
    p <- p + ggplot2::geom_vline(xintercept = x$fema$b)
  } else {
    p <- p + ggplot2::geom_vline(xintercept = ref)
  }
  
  p <- p + ggplot2::scale_x_continuous(
    name   = xlab,
    limits = range(pretty(range(xlim))),
    breaks = pretty(range(xlim)))
  if (x$method %in% c("P-FPV", "D-FAV", "D-FIV")){
    p <- p + ggplot2::scale_y_reverse(name = ylab, breaks = 1/yax, labels = yax, limits = rev(range(pretty(ylim))))
  } else if (x$method %in% c("E-UW", "E-FIV")){
    p <- p + ggplot2::scale_y_reverse(name = ylab, limits = rev(range(pretty(ylim))), breaks = pretty(range(ylim)))
  } else {
    p <- p + ggplot2::scale_y_continuous(name = ylab, limits = range(pretty(ylim)), breaks = pretty(range(ylim)))
  }
  
  return(p)
}
