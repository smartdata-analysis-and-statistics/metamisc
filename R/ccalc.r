#' Calculate the concordance statistic
#'
#' The function calculates (transformed versions of) the concordance (c-)
#' statistic with the corresponding sampling variance. 
#' 
#' @param cstat vector to specify the estimated c-statistics.
#' @param cstat.se Optional vector to specify the corresponding standard errors.
#' @param cstat.cilb Optional vector to specify the lower limits of the confidence interval.
#' @param cstat.ciub Optional vector to specify the upper limits of the confidence interval.
#' @param cstat.cilv Optional vector to specify the levels of aformentioned confidence interval limits. 
#' (default: 0.95, which corresponds to the 95\% confidence interval).
#' @param sd.LP Optional vector to specify the standard deviations of the linear predictor (prognostic index).
#' @param N Optional vector to specify the sample/group sizes.
#' @param O Optional vector to specify the total number of observed events.
#' @param Po Optional vector to specify the observed event probabilities.
#' @param data Optional data frame containing the variables given to the arguments above.
#' @param slab Optional vector with labels for the studies.
#' @param subset Optional vector indicating the subset of studies that should be used. This can be a logical vector or a numeric vector indicating the indices of the studies to include.
#' @param g a quoted string that is the function to transform estimates of the c-statistic; see the details below.
#' @param level Optional numeric to specify the level for the confidence interval, default \code{0.95}.
#' @param approx.se.method integer specifying which method should be used for estimating the standard error of the
#' c-statistic (Newcombe, 2006). So far, only method \code{2} and method \code{4} (default) have been implemented.
#' @param \ldots Additional arguments.
#' 
#' @details 
#' The c-statistic is a measure of discrimination, and indicates the ability of a prediction model to 
#' distinguish between patients developing and not developing the outcome. The c-statistic typically ranges 
#' from 0.5 (no discriminative ability) to 1 (perfect discriminative ability). 
#' 
#' By default, the function \code{ccalc} will derive the c-statistic of each study, together with 
#' the corresponding standard error and 95\% confidence interval. However, it is also possible to calculate transformed 
#' versions of the c-statistic. Appropriate standard errors are then derived using the Delta method. 
#' For instance, the logit transformation can be applied by specifying \code{g="log(cstat/(1-cstat))"}. 
#' 
#' \subsection{Restoring the c-statistic}{
#' For studies where the c-statistic is missing, it is estimated from the standard deviation of the linear predictor 
#' (\code{theta.source="std.dev(LP)"}). The corresponding method is described by White et al. (2015).
#' }
#' 
#' \subsection{Restoring the standard error of the c-statistic}{
#' When missing, the standard error of the c-statistic can be estimated from the confidence interval. Alternatively, 
#' the standard error can be approximated from a combination of the reported c-statistic, the total sample size and 
#' the total number of events (Newcombe, 2006). This can be achieved by adopting (a modification of) the method 
#' proposed by Hanley and McNeil, as specified in \code{approx.se.method}.
#' }
#' 
#' @references 
#' Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review and meta-analysis of prediction model performance. BMJ. 2017;356:i6460. 
#' 
#' Debray TPA, Damen JAAG, Riley R, Snell KIE, Reitsma JB, Hooft L, et al. A framework for meta-analysis of  prediction model studies with binary and time-to-event outcomes. Stat Methods Med Res. 2018; In press. 
#' 
#' Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) 
#' curve. \emph{Radiology}. 1982; 143(1):29--36.
#' 
#' Newcombe RG. Confidence intervals for an effect size measure based on the Mann-Whitney statistic. 
#' Part 2: asymptotic methods and evaluation. \emph{Stat Med}. 2006; 25(4):559--73.
#' 
#' Snell KI, Ensor J, Debray TP, Moons KG, Riley RD. Meta-analysis of prediction model performance across 
#' multiple studies: Which scale helps ensure between-study normality for the C -statistic and calibration measures? 
#' \emph{Statistical Methods in Medical Research}. 2017. 
#' 
#' White IR, Rapsomaniki E, the Emerging Risk Factors Collaboration. Covariate-adjusted measures of discrimination 
#' for survival data. \emph{Biom J}. 2015;57(4):592--613. 
#' 
#' 
#' @return An object of class c("mm_perf","data.frame") with the following columns:
#' \describe{
##'  \item{"theta"}{The (transformed) c-statistics. }
##'  \item{"theta.se"}{Standard errors of the (transformed) c-statistics.}
##'  \item{"theta.cilb"}{Lower confidence interval of the (transformed) c-statistics. The level is specified in
##'  \code{level}. Intervals are calculated on the same scale as \code{theta} by assuming a Normal distribution.}
##'  \item{"theta.ciub"}{Upper confidence interval of the (transformed) c-statistics. The level is specified in
##'  \code{level}. Intervals are calculated on the same scale as \code{theta} by assuming a Normal distribution.}
##'  \item{"theta.source"}{Method used for calculating the (transformed) c-statistic.}
##'  \item{"theta.se.source"}{Method used for calculating the standard error of the (transformed) c-statistic.}
##' }
#' 
#' @examples 
#' ######### Validation of prediction models with a binary outcome #########
#' data(EuroSCORE)
#' 
#' # Calculate the c-statistic and its standard error
#' est1 <- ccalc(cstat = c.index, cstat.se = se.c.index, cstat.cilb = c.index.95CIl, 
#'               cstat.ciub = c.index.95CIu, N = n, O = n.events, data = EuroSCORE, slab = Study)
#' est1
#'   
#' # Calculate the logit c-statistic and its standard error
#' est2 <- ccalc(cstat = c.index, cstat.se = se.c.index, cstat.cilb = c.index.95CIl, 
#'               cstat.ciub = c.index.95CIu, N = n, O = n.events, data = EuroSCORE, slab = Study, 
#'               g = "log(cstat/(1-cstat))")
#' est2
#'       
#' # Display the results of all studies in a forest plot
#' plot(est1)
#'                                                             
#' @keywords meta-analysis discrimination concordance statistic performance
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
#' 
ccalc <- function(cstat, cstat.se, cstat.cilb, cstat.ciub, cstat.cilv, sd.LP, N, O, Po, data, slab, subset,
                  g = NULL, level = 0.95, approx.se.method = 4, ...) {
  
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
  slab          <- eval(mf.slab,   data, enclos=sys.frame(sys.parent()))
  mf.subset     <- mf[[match("subset", names(mf))]]
  subset        <- eval(mf.subset, data, enclos=sys.frame(sys.parent()))
  mf.cstat      <- mf[[match("cstat", names(mf))]]
  cstat         <- eval(mf.cstat, data, enclos=sys.frame(sys.parent()))
  mf.cstat.se   <- mf[[match("cstat.se", names(mf))]]
  cstat.se      <- eval(mf.cstat.se, data, enclos=sys.frame(sys.parent()))
  mf.cstat.cilb <- mf[[match("cstat.cilb", names(mf))]]
  cstat.cilb    <- eval(mf.cstat.cilb, data, enclos=sys.frame(sys.parent()))
  mf.cstat.ciub <- mf[[match("cstat.ciub", names(mf))]]
  cstat.ciub    <- eval(mf.cstat.ciub, data, enclos=sys.frame(sys.parent()))
  mf.cstat.cilv <- mf[[match("cstat.cilv", names(mf))]]
  cstat.cilv    <- eval(mf.cstat.cilv, data, enclos=sys.frame(sys.parent()))
  mf.sd.LP      <- mf[[match("sd.LP", names(mf))]]
  sd.LP         <- eval(mf.sd.LP, data, enclos=sys.frame(sys.parent()))
  mf.N          <- mf[[match("N", names(mf))]]
  N             <- eval(mf.N, data, enclos=sys.frame(sys.parent()))
  mf.O          <- mf[[match("O", names(mf))]]
  O             <- eval(mf.O, data, enclos=sys.frame(sys.parent()))
  mf.Po         <- mf[[match("Po", names(mf))]]
  Po            <- eval(mf.Po, data, enclos=sys.frame(sys.parent()))
  
  
  
  #######################################################################################
  # Count number of studies
  #######################################################################################
  k <- 0
  
  if (!no.data) {
    k <- dim(data)[1]
  } else if (!is.null(cstat)) {
    k <- length(cstat)
  } else if (!is.null(cstat.se)) {
    k <- length(cstat.se)
  } else if (!is.null(cstat.cilb)) {
    k <- length(cstat.cilb)
  } else if (!is.null(cstat.ciub)) {
    k <- length(cstat.ciub)
  } else if (!is.null(sd.LP)) {
    k <- length(sd.LP)
  } else if (!is.null(N)) {
    k <- length(N)
  } else if (!is.null(O)) {
    k <- length(O)
  } else if (!is.null(Po)) {
    k <- length(Po)
  }
  
  if (k<1) stop("No data provided!")
  
  if(is.null(cstat)) {
    cstat <- rep(NA, times=k)
  }
  
  #######################################################################################
  # Assign study labels
  # taken from escalc
  #######################################################################################
  if (!is.null(slab)) {
    
    if (!is.null(subset))
      slab <- slab[subset]
    
    if (anyNA(slab))
      stop("NAs in study labels.")
    
    if (inherits(slab, "factor")) {
      slab <- as.character(slab)
    }
    
    ### check if study labels are unique; if not, make them unique
    
    if (anyDuplicated(slab))
      slab <- make.unique(slab)
    
    if (length(slab) != k)
      stop("Study labels not of same length as data.")
    
    ### add slab attribute to the cstat vector
    attr(cstat, "slab") <- slab
  }
  
  ### if a subset of studies is specified (note: subsetting of other parts already done above, so yi/vi/ni.u/slab are already subsetted)
  if (!is.null(subset)) {
    if (!no.data)
      data <- data[subset,,drop=FALSE]
  }
  
  #######################################################################################
  # Prepare data
  #######################################################################################
  if (is.null(O))  O <- rep(NA, length=k)
  if (is.null(Po)) Po <- rep(NA, length=k)
  if (is.null(N))  N <- rep(NA, length=k)
  if (is.null(cstat.cilb)) cstat.cilb <- rep(NA, length=k)
  if (is.null(cstat.ciub)) cstat.ciub <- rep(NA, length=k)
  if (is.null(cstat.cilv)) cstat.cilv <- rep(0.95, length=k)
  if (is.null(cstat.se)) cstat.se <- array(NA, dim=k)
  if (is.null(sd.LP)) sd.LP <- rep(NA, k)
  
  if (sum(cstat.cilv>1 | cstat.cilv<0) > 0) {
    stop("Invalid level(s) specified for 'cstat.cilv'!")
  }
  if (is.na(level) | level<0 | level>1) {
    stop("Invalid value for 'level'!")
  }
  if (!approx.se.method %in% c(2,4)) {
    stop("Invalid method for restoring the SE of the c-statistic!")
  }
  
  # Calculate O and N from other information if possible
  O <- ifelse(is.na(O), Po*N, O)
  N <- ifelse(is.na(N), O/Po, N)
  
  theta.cil <- theta.ciu <- rep(NA, k)
  
  # Restore c-statistic
  te.method <- c("c-statistic", "std.dev(LP)")
  te.orig   <- calculate.cstat.theta(cstat=cstat, g=g)
  te.white  <- calculate.cstat.sdPI(sdPI=sd.LP, g=g)
  
  te.dat <- cbind(te.orig, te.white)
  
  # For each study, find the first colum without missing
  myfun = function(dat) { which.min(is.na(dat)) }
  
  sel.cstat <- apply(te.dat, 1, myfun)
  theta <- te.dat[cbind(seq_along(sel.cstat), sel.cstat)]                            
  theta.source <-  te.method[sel.cstat]
  
  # Directly transform the provided confidence limits for studies where the reported level is equal to the requested level
  if (is.null(g)) {
    theta.cil[cstat.cilv==level] <- cstat.cilb[cstat.cilv==level]
    theta.ciu[cstat.cilv==level] <- cstat.ciub[cstat.cilv==level]
  } else {
    theta.cil[cstat.cilv==level] <- eval(parse(text=g), list(cstat = cstat.cilb[cstat.cilv==level]))
    theta.ciu[cstat.cilv==level] <- eval(parse(text=g), list(cstat = cstat.ciub[cstat.cilv==level]))
  } 
  
  # Calculate all the possible variations of var(theta)
  tv.method  <- c("Standard Error", "Confidence Interval", paste("Newcombe (Method ", approx.se.method, ")", sep=""), paste("Newcombe (Method ", approx.se.method, ")", sep=""))
  tv.se      <- restore.c.var.se(cstat=cstat, c.se=cstat.se, g=g) # Derived from standard error
  tv.ci      <- restore.c.var.ci(cil=cstat.cilb, ciu=cstat.ciub, level=cstat.cilv, g=g) # Derived from 95% confidence interval
  tv.hanley  <- restore.c.var.hanley(cstat=cstat, N.subjects=N, N.events=O, restore.method=approx.se.method, g=g)
  tv.hanley2 <- restore.c.var.hanley2(sd.LP=sd.LP, N.subjects=N, N.events=O, restore.method=approx.se.method, g=g)
  
  # Save all estimated variances. The order of the columns indicates the priority             
  dat <-cbind(tv.se, tv.ci, tv.hanley, tv.hanley2)  
  
  sel.var <- apply(dat, 1, myfun)
  theta.var <- dat[cbind(seq_along(sel.var), sel.var)]                            
  theta.var.source <-  tv.method[sel.var]
  
  # Calculate the desired confidence intervals
  theta.cil[is.na(theta.cil)] <- (theta+qnorm((1-level)/2)*sqrt(theta.var))[is.na(theta.cil)]
  theta.ciu[is.na(theta.ciu)] <- (theta+qnorm((1+level)/2)*sqrt(theta.var))[is.na(theta.ciu)]
  
  
  # Store results, and method for calculating SE
  ds <- data.frame(theta=theta, theta.se=sqrt(theta.var), theta.cilb=theta.cil, theta.ciub=theta.ciu, 
                   theta.source=theta.source, theta.se.source=theta.var.source)
  
  if(is.null(slab) & !no.data) {
    slab <- rownames(data)
    rownames(ds) <- slab
  } else if (!is.null(slab)) {
    slab <- make.unique(as.character(slab))
    rownames(ds) <- slab
  }
  
  # Add some attributes specifying the nature of the (untransformed) estimatess
  attr(ds, 'estimand') <- "c-statistic"
  attr(ds, 'theta_scale') <- g
  attr(ds, 'plot_refline') <- 0.5
  attr(ds, 'plot_lim') <- c(0,1)
  
  class(ds) <- c("mm_perf", class(ds))
  
  return(ds)
}