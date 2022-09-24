#' Calculate the total O:E ratio
#'
#' This function calculates (transformed versions of) the ratio of total number of observed versus expected events with the 
#' corresponding sampling variance. 
#' 
#' @param OE vector with the estimated ratio of total observed versus total expected events
#' @param OE.se Optional vector with the standard errors of the estimated O:E ratios. 
#' @param OE.cilb Optional vector to specify the lower limits of the confidence interval for \code{OE}.
#' @param OE.ciub Optional vector to specify the upper limits of the confidence interval for \code{OE}.
#' @param OE.cilv Optional vector to specify the levels of aformentioned confidence interval limits. 
#' (default: 0.95, which corresponds to the 95\% confidence interval).
#' @param EO Optional vector with the estimated ratio of total expected versus total observed events
#' @param EO.se Optional vector with the standard errors of the estimated E:O ratios
#' @param citl Optional vector with the estimated calibration-in-the-large statistics
#' @param citl.se Optional vector with the standard error of the calibration-in-the-large statistics
#' @param N Optional vector to specify the sample/group sizes.
#' @param O Optional vector to specify the total number of observed events.
#' @param E Optional vector to specify the total number of expected events
#' @param Po Optional vector to specify the (cumulative) observed event probabilities.
#' @param Po.se Optional vector with the standard errors of \code{Po}. For time-to-event data, 
#' these could also be the SE of the observed survival probabilities (e.g. as obtained from Kaplan-Meier analysis)
#' @param Pe Optional vector to specify the (cumulative) expected event probabilites
#' (if specified, during time \code{t.val})
#' @param data Optional data frame containing the variables given to the arguments above.
#' @param slab Optional vector with labels for the studies.
#' @param add a non-negative number indicating the amount to add to zero counts. See `Details'
#' @param g a quoted string that is the function to transform estimates of the total O:E ratio; see the details below.
#' @param level level for confidence interval, default \code{0.95}.
#' @param \ldots Additional arguments.
#' 
#' @references 
#' Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review and meta-analysis of prediction model performance. \emph{BMJ}. 2017;356:i6460. 
#' 
#' Debray TPA, Damen JAAG, Riley R, Snell KIE, Reitsma JB, Hooft L, et al. A framework for meta-analysis of prediction model studies with binary and time-to-event outcomes. \emph{Stat Methods Med Res}. 2019 Sep;28(9):2768--86.
#'  
#' Snell KI, Ensor J, Debray TP, Moons KG, Riley RD. Meta-analysis of prediction model performance across 
#' multiple studies: Which scale helps ensure between-study normality for the C -statistic and calibration measures? 
#' \emph{Stat Methods Med Res}. 2017. 
#' 
#' 
#' @return An object of class c("mm_perf","data.frame") with the following columns:
#' \describe{
##'  \item{"theta"}{The (transformed) O:E ratio. }
##'  \item{"theta.se"}{Standard errors of the (transformed) O:E ratio.}
##'  \item{"theta.cilb"}{Lower confidence interval of the (transformed) O:E ratios. The level is specified in
##'  \code{level}. Intervals are calculated on the same scale as \code{theta} by assuming a Normal distribution.}
##'  \item{"theta.ciub"}{Upper confidence interval of the (transformed) c-statistics. The level is specified in
##'  \code{level}. Intervals are calculated on the same scale as \code{theta} by assuming a Normal distribution.}
##'  \item{"theta.source"}{Method used for calculating the (transformed) O:E ratio.}
##'  \item{"theta.se.source"}{Method used for calculating the standard error of the (transformed) O:E ratio.}
##' }
#' 
#' @examples 
#' ######### Validation of prediction models with a binary outcome #########
#' data(EuroSCORE)
#' 
#' # Calculate the total O:E ratio and its standard error
#' est1 <- oecalc(O = n.events, E = e.events, N = n, data = EuroSCORE, slab = Study)
#' est1
#' 
#' # Calculate the log of the total O:E ratio and its standard error
#' est2 <- oecalc(O = n.events, E = e.events, N = n, data = EuroSCORE, slab = Study, g = "log(OE)")
#' est2
#' 
#' # Display the results of all studies in a forest plot
#' plot(est1)
#' 
#' @keywords meta-analysis calibration performance
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
oecalc <- function(OE, OE.se, OE.cilb, OE.ciub, OE.cilv, EO, EO.se, citl, citl.se, N, O, E, Po, Po.se, Pe, 
                   data, slab, add=1/2, g=NULL, level=0.95, ...) {
  
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
  mf.OE         <- mf[[match("OE", names(mf))]]
  OE            <- eval(mf.OE, data, enclos=sys.frame(sys.parent()))
  mf.OE.se      <- mf[[match("OE.se", names(mf))]]
  OE.se         <- eval(mf.OE.se, data, enclos=sys.frame(sys.parent()))
  mf.OE.cilb    <- mf[[match("OE.cilb", names(mf))]]
  OE.cilb       <- eval(mf.OE.cilb, data, enclos=sys.frame(sys.parent()))
  mf.OE.ciub    <- mf[[match("OE.ciub", names(mf))]]
  OE.ciub       <- eval(mf.OE.ciub, data, enclos=sys.frame(sys.parent()))
  mf.OE.cilv    <- mf[[match("OE.cilv", names(mf))]]
  OE.cilv       <- eval(mf.OE.cilv, data, enclos=sys.frame(sys.parent()))
  mf.EO         <- mf[[match("EO", names(mf))]]
  EO            <- eval(mf.EO, data, enclos=sys.frame(sys.parent()))
  mf.EO.se      <- mf[[match("EO.se", names(mf))]]
  EO.se         <- eval(mf.EO.se, data, enclos=sys.frame(sys.parent()))
  mf.citl       <- mf[[match("citl", names(mf))]]
  citl          <- eval(mf.citl, data, enclos=sys.frame(sys.parent()))
  mf.citl.se    <- mf[[match("citl.se", names(mf))]]
  citl.se       <- eval(mf.citl.se, data, enclos=sys.frame(sys.parent()))
  mf.N          <- mf[[match("N", names(mf))]]
  N             <- eval(mf.N, data, enclos=sys.frame(sys.parent()))
  mf.O          <- mf[[match("O", names(mf))]]
  O             <- eval(mf.O, data, enclos=sys.frame(sys.parent()))
  mf.E          <- mf[[match("E", names(mf))]]
  E             <- eval(mf.E, data, enclos=sys.frame(sys.parent()))
  mf.Po         <- mf[[match("Po", names(mf))]]
  Po            <- eval(mf.Po, data, enclos=sys.frame(sys.parent()))
  mf.Pe         <- mf[[match("Pe", names(mf))]]
  Pe            <- eval(mf.Pe, data, enclos=sys.frame(sys.parent()))
  mf.Po.se      <- mf[[match("Po.se", names(mf))]]
  Po.se         <- eval(mf.Po.se, data, enclos=sys.frame(sys.parent()))
  
  
  #######################################################################################
  # Count number of studies
  #######################################################################################
  k <- 0
  
  if (!no.data) {
    k <- dim(data)[1]
  } else if (!is.null(OE)) {
    k <- length(OE)
  } else if (!is.null(OE.se)) {
    k <- length(OE.se)
  } else if (!is.null(EO)) {
    k <- length(EO)
  } else if (!is.null(EO.se)) {
    k <- length(EO.se)
  } else if (!is.null(OE.cilb)) {
    k <- length(OE.cilb)
  } else if (!is.null(OE.ciub)) {
    k <- length(OE.ciub)
  } else if (!is.null(OE.cilv)) {
    k <- length(OE.cilv)
  } else if (!is.null(citl)) {
    k <- length(citl)
  } else if (!is.null(citl.se)) {
    k <- length(citl.se)
  } else if (!is.null(N)) {
    k <- length(N)
  }  else if (!is.null(O)) {
    k <- length(O)
  } else if (!is.null(E)) {
    k <- length(E)
  } else if (!is.null(Po)) {
    k <- length(Po)
  } else if (!is.null(Po.se)) {
    k <- length(Po.se)
  } else if (!is.null(Pe)) {
    k <- length(Pe)
  }

  if (k<1) stop("No data provided!")
  
  if(is.null(OE))  OE <- rep(NA, times=k)
  if(is.null(OE.se)) OE.se <- rep(NA, times=k)
  if(is.null(EO))  EO <- rep(NA, times=k)
  if(is.null(EO.se)) EO.se <- rep(NA, times=k)
  if(is.null(OE.cilb)) OE.cilb <- rep(NA, times=k)
  if(is.null(OE.ciub)) OE.ciub <- rep(NA, times=k)
  if(is.null(OE.cilv)) OE.cilv <- rep(0.95, times=k) # Assume 95% CI by default 
  if(is.null(O)) O <- rep(NA, times=k)
  if(is.null(E)) E <- rep(NA, times=k)
  if(is.null(N)) N <- rep(NA, times=k)
  if(is.null(Po)) Po <- rep(NA, times=k)
  if(is.null(Pe)) Pe <- rep(NA, times=k)
  if(is.null(Po.se)) Po.se <- rep(NA, times=k)
  
  
  #######################################################################################
  # Assign study labels
  # taken from escalc
  #######################################################################################
  if (!is.null(slab)) {
    
    if (anyNA(slab))
      stop("NAs in study labels.")
    
    if (class(slab)=="factor") {
      slab <- as.character(slab)
    }
    
    ### check if study labels are unique; if not, make them unique
    if (anyDuplicated(slab))
      slab <- make.unique(slab)
    
    if (length(slab) != k)
      stop("Study labels not of same length as data.")
    
    ### add slab attribute to the cstat vector
    attr(OE, "slab") <- slab
  }
  

  #######################################################################################
  # Derive the OE ratio and its error variance
  # The order defines the preference for the final result
  #######################################################################################
  results <- list()
  results[[1]] <- data.frame(est=resoe.OE.se(OE=OE, OE.se=OE.se, g=g), method="OE and SE(OE)")
  results[[2]] <- data.frame(est=resoe.OE.ci(OE=OE, OE.cilb=OE.cilb, OE.ciub=OE.ciub, OE.cilv=OE.cilv, g=g), method="OE and CI(OE)")                                   #                                model=pars.default$model.oe)
  results[[3]] <- data.frame(est=resoe.O.E.N(O=O, E=E, N=N, correction = add, g=g), method="O, E and N")
  results[[4]] <- data.frame(est=resoe.O.Pe.N(O=O, Pe=Pe, N=N, correction = add, g=g), method="O, Pe and N")
  results[[5]] <- data.frame(est=resoe.E.Po.N(E=E, Po=Po, N=N, correction = add, g=g), method="E, Po and N")
  results[[6]] <- data.frame(est=resoe.Po.Pe.N(Po=Po, Pe=Pe, N=N, correction = add, g=g), method="Po, Pe and N")
  results[[7]] <- data.frame(est=resoe.EO.se(EO=EO, EO.se=EO.se, g=g), method="EO and SE(EO)")
  results[[8]] <- data.frame(est=resoe.O.Po.E(O=O, Po=Po, E=E, correction = add, g=g), method="O, E and Po") 
  results[[9]] <- data.frame(est=resoe.O.Pe.E(O=O, Pe=Pe, E=E, correction = add, g=g), method="O, E and Pe") 
  results[[10]] <- data.frame(est=resoe.O.E(O=O, E=E, correction = add, g=g), method="O and E")
  results[[11]] <- data.frame(est=resoe.Po.Pe.sePo(Po=Po, Pe=Pe, Po.se=Po.se, g=g), method = "Po, Pe and SE(Po)")
  results[[12]] <- data.frame(est=resoe.Po.Pe(Po=Po, Pe=Pe, g=g), method = "Po and Pe")


  #t.citl   <- restore.oe.citl(citl=citl, citl.se=citl.se, O=O, Po=Po, N=N, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, 
  #                            model=pars.default$model.oe) 
  
  # Select appropriate estimate for 'theta' and record its source
  dat.est <- dat.se <- dat.method <-NULL
  for (i in 1:length(results)) {
    dat.est <- cbind(dat.est, (results[[i]])[,1]) 
    dat.se <- cbind(dat.se, sqrt(results[[i]][,2])) #take square root of error variance
    dat.method <- cbind(dat.method, as.character(results[[i]][,3]))
  }

  myfun = function(dat) { which.min(is.na(dat)) }
  sel.theta1 <- apply(dat.est, 1, myfun)
  sel.theta2 <- apply(dat.se, 1, myfun)
  sel.theta <- apply(cbind(sel.theta1, sel.theta2), 1, max) # only take the estimate for which we have an SE available
  
  theta <- dat.est[cbind(seq_along(sel.theta), sel.theta)] # Preferred estimate for theta 
  theta.se <- dat.se[cbind(seq_along(sel.theta), sel.theta)]# Preferred estimate for SE(theta)
  theta.source <-  dat.method[cbind(seq_along(sel.theta), sel.theta)] # Method used for estimating theta and its SE
  
  #######################################################################################
  # Derive confindence intervals
  #######################################################################################
  theta.cil <- theta.ciu <- rep(NA, k)
  
  # Directly transform the provided confidence limits for studies where the reported level is equal to the requested level
  if (is.null(g)) {
    theta.cil[OE.cilv==level] <- OE.cilb[OE.cilv==level]
    theta.ciu[OE.cilv==level] <- OE.ciub[OE.cilv==level]
  } else {
    theta.cil[OE.cilv==level] <- eval(parse(text=g), list(OE = OE.cilb[OE.cilv==level]))
    theta.ciu[OE.cilv==level] <- eval(parse(text=g), list(OE = OE.ciub[OE.cilv==level]))
  } 
  
  # Calculate the desired confidence intervals
  theta.cil[is.na(theta.cil)] <- (theta+qnorm((1-level)/2)*theta.se)[is.na(theta.cil)]
  theta.ciu[is.na(theta.ciu)] <- (theta+qnorm((1+level)/2)*theta.se)[is.na(theta.ciu)]
  
  #######################################################################################
  # Attempt to restore O, E and N
  #######################################################################################
  O[is.na(O)] <- (Po*N)[is.na(O)]
  E[is.na(E)] <- (Pe*N)[is.na(E)]
  N[is.na(N)] <- (O/Po)[is.na(N)]
  N[is.na(N)] <- (E/Pe)[is.na(N)]
  
  
  #######################################################################################
  # Sore results
  #######################################################################################
  ds <- data.frame(theta=theta, theta.se=theta.se, 
                   theta.cilb=theta.cil, theta.ciub=theta.ciu, 
                   theta.source=theta.source, O=O, E=E, N=N)
  
  if(is.null(slab) & !no.data) {
    slab <- rownames(data)
    rownames(ds) <- slab
  } else if (!is.null(slab)) {
    slab <- make.unique(as.character(slab))
    rownames(ds) <- slab
  }
  
  # Add some attributes specifying the nature of the data
  attr(ds, 'estimand') <- "Total O:E ratio"
  attr(ds, 'theta_scale') <- g
  attr(ds, 'plot_refline') <- 1
  attr(ds, 'plot_lim') <- c(0,NA)
  
  class(ds) <- c("mm_perf", class(ds))
  
  return(ds)
}


resoe.O.Pe.N <- function(O, Pe, N, correction, g=NULL) {
  return(resoe.O.E.N(O=O, E=Pe*N, N=N, correction=correction, g=g))
}

resoe.E.Po.N <- function(E, Po, N, correction, g=NULL) {
  return(resoe.O.E.N(O=Po*N, E=E, N=N, correction=correction, g=g))
}

resoe.Po.Pe.N <- function(Po, Pe, N, correction, g=NULL) {
  return(resoe.O.E.N(O=Po*N, E=Pe*N, N=N, correction=correction, g=g))
}

resoe.O.Po.E <- function(O, Po, E, correction, g=NULL) {
  return(resoe.O.E.N(O=O, E=E, N=O/Po, correction=correction, g=g))
}

resoe.O.Pe.E <- function(O, Pe, E, correction, g=NULL) {
  return(resoe.O.E.N(O=O, E=E, N=E/Pe, correction=correction, g=g))
}


# Calculate OE and its error variance from O, E and N
resoe.O.E.N <- function(O, E, N, correction, g=NULL) {
  numstudies <- length(O)
  dat_OE <- dat_g_OE <- matrix(NA, nrow = numstudies, ncol = 2)
  colnames(dat_OE) <- c("OE", "var_OE")
  colnames(dat_g_OE) <- c("est", "var_est")
  
  cc <- which(E == 0)
  E[cc] <- E[cc] + correction
  O[cc] <- O[cc] + correction
  N[cc] <- N[cc] + correction
  
  dat_OE[,1] <- O/E 
  dat_OE[,2] <- ((O*(1-O/N))/(E**2)) # Error variance
  
  if (!is.null(g)) {
    # Apply delta method to the transformation
    gd <- deriv(parse(text = g), "OE") # Take the derivative of g
    
    dat <- data.frame(OE = dat_OE[,1])
    g_OE <- eval(gd, envir = dat)
    dat_g_OE[,1] <-  g_OE # Extract the transformed OE ratio
    dat_g_OE[,2] <- (attr(g_OE, "grad")**2) * dat_OE[,2] # Derive the variance of the transformed OE ratio
    return(dat_g_OE)
  }
  return(dat_OE)
}

resoe.Po.Pe.sePo <- function(Po, Pe, Po.se, g = NULL) {
  numstudies <- length(Po)
  dat_OE <- dat_g_OE <- matrix(NA, nrow = numstudies, ncol = 2)
  colnames(dat_OE) <- c("OE", "var_OE")
  colnames(dat_g_OE) <- c("est", "var_est")
  
  dat_OE[,1] <- Po/Pe
  dat_OE[,2] <- Po.se**2/(Pe**2)
  
  if (!is.null(g)) {
    # Apply delta method to the transformation
    gd <- deriv(parse(text = g), "OE") # Take the derivative of g
    
    dat <- data.frame(OE = dat_OE[,1])
    g_OE <- eval(gd, envir = dat)
    dat_g_OE[,1] <-  g_OE # Extract the transformed OE ratio
    dat_g_OE[,2] <- (attr(g_OE, "grad")**2) * dat_OE[,2] # Derive the variance of the transformed OE ratio
    return(dat_g_OE)
  }
  return(dat_OE)
}

resoe.EO.se  <- function(EO, EO.se, g = NULL) {
  numstudies <- length(EO)
  dat_OE <- dat_g_OE <- matrix(NA, nrow = numstudies, ncol = 2)
  colnames(dat_OE) <- c("OE", "var_OE")
  colnames(dat_g_OE) <- c("est", "var_est")
  
  dat_OE[,1] <- 1/EO
  dat_OE[,2] <- (EO.se**2)/(EO**4)
 
  if (!is.null(g)) {
    # Apply delta method to the transformation
    gd <- deriv(parse(text = g), "OE") # Take the derivative of g
    
    dat <- data.frame(OE = dat_OE[,1])
    g_OE <- eval(gd, envir = dat)
    dat_g_OE[,1] <-  g_OE # Extract the transformed OE ratio
    dat_g_OE[,2] <- (attr(g_OE, "grad")**2) * dat_OE[,2] # Derive the variance of the transformed OE ratio
    return(dat_g_OE)
  }
  return(dat_OE)
}
  
resoe.OE.se <- function(OE, OE.se, g = NULL) {
  numstudies <- length(OE)
  dat_OE <- dat_g_OE <- matrix(NA, nrow = numstudies, ncol = 2)
  colnames(dat_OE) <- c("OE", "var_OE")
  colnames(dat_g_OE) <- c("est", "var_est")
  
  dat_OE[,1] <- OE
  dat_OE[,2] <- OE.se**2
  
  if (!is.null(g)) {
    # Apply delta method to the transformation
    gd <- deriv(parse(text = g), "OE") # Take the derivative of g
    
    dat <- data.frame(OE = dat_OE[,1])
    g_OE <- eval(gd, envir = dat)
    dat_g_OE[,1] <-  g_OE # Extract the transformed OE ratio
    dat_g_OE[,2] <- (attr(g_OE, "grad")**2) * dat_OE[,2] # Derive the variance of the transformed OE ratio
    return(dat_g_OE)
  }
  return(dat_OE)
}

resoe.OE.ci <- function(OE, OE.cilb, OE.ciub, OE.cilv, g=NULL) {
  numstudies <- length(OE)
  dat_OE <- dat_g_OE <- matrix(NA, nrow = numstudies, ncol = 2)
  colnames(dat_OE) <- c("OE", "var_OE")
  colnames(dat_g_OE) <- c("est", "var_est")
  
  dat_OE[,1] <- OE
  dat_OE[,2] <- ((OE.ciub - OE.cilb)/(2*qnorm(0.5 + OE.cilv/2)))**2 #Derive the error variance from 95% CI
  
  if (!is.null(g)) {
    # Apply delta method to the transformation
    gd <- deriv(parse(text = g), "OE") # Take the derivative of g
    
    dat <- data.frame(OE = dat_OE[,1])
    g_OE <- eval(gd, envir = dat)
    dat_g_OE[,1] <-  g_OE # Extract the transformed OE ratio
    dat_g_OE[,2] <- (attr(g_OE, "grad")**2) * dat_OE[,2] # Derive the variance of the transformed OE ratio
    return(dat_g_OE)
  }
  return(dat_OE)
}

resoe.O.E <- function(O, E, correction, g = NULL) {
  numstudies <- length(O)
  dat_OE <- dat_g_OE <- matrix(NA, nrow = numstudies, ncol = 2)
  colnames(dat_OE) <- c("OE", "var_OE")
  colnames(dat_g_OE) <- c("est", "var_est")
  
  cc <- which(E == 0)
  E[cc] <- E[cc] + correction
  O[cc] <- O[cc] + correction
  
  dat_OE[,1] <- O/E
  dat_OE[,2] <- (O/(E**2))
  
  if (!is.null(g)) {
    # Apply delta method to the transformation
    gd <- deriv(parse(text = g), "OE") # Take the derivative of g
    
    dat <- data.frame(OE = dat_OE[,1])
    g_OE <- eval(gd, envir = dat)
    dat_g_OE[,1] <-  g_OE # Extract the transformed OE ratio
    dat_g_OE[,2] <- (attr(g_OE, "grad")**2) * dat_OE[,2] # Derive the variance of the transformed OE ratio
    return(dat_g_OE)
  }
  return(dat_OE)
}

resoe.citl <- function(citl, citl.se, Po, O, N, correction, g = NULL) {
  numstudies <- length(O)
  dat_citl <- matrix(NA, nrow = numstudies, ncol = 2)
  colnames(dat_citl) <- c("citl", "var_citl")
  
  warning("Implementation not finalized!")
  return(dat_citl)
}

resoe.Po.Pe <- function (Po, Pe, g = NULL) {
  numstudies <- length(Po)
  dat_OE <- dat_g_OE <- matrix(NA, nrow = numstudies, ncol = 2)
  colnames(dat_OE) <- c("OE", "var_OE")
  colnames(dat_g_OE) <- c("est", "var_est")
  
  dat_OE[,1] <- Po/Pe
  dat_OE[,2] <- NA # SE cannot be estimated from Po and Pe alone
  
  if (!is.null(g)) {
    # Apply delta method to the transformation
    gd <- deriv(parse(text = g), "OE") # Take the derivative of g
    
    dat <- data.frame(OE = dat_OE[,1])
    g_OE <- eval(gd, envir = dat)
    dat_g_OE[,1] <-  g_OE # Extract the transformed OE ratio
    dat_g_OE[,2] <- NA  # SE cannot be estimated from Po and Pe alone
    return(dat_g_OE)
  }
  return(dat_OE)
}
