

# Adapted from R package 'car' to avoid package loading issues in R-forge
deltaMethod <- function (object, g, vcov., func = g, constants, level = 0.95, ...) {
  if (!is.character(g)) 
    stop("The argument 'g' must be a character string")
  
  para <- object         
  g <- parse(text = g)
  q <- length(para)
  for (i in 1:q) {
    assign(names(para)[i], para[i])
  }
  if(!missing(constants)){
    for (i in seq_along(constants)) assign(names(constants[i]), constants[[i]])}
  est <- eval(g)
  names(est) <- NULL
  gd <- rep(0, q)
  for (i in 1:q) {
    gd[i] <- eval(D(g, names(para)[i]))
  }
  se.est <- as.vector(sqrt(t(gd) %*% vcov. %*% gd))
  result <- data.frame(Estimate = est, SE = se.est, row.names = c(func))
  result
}

# Adapted from matrixcalc
hadamard.prod <- function (x, y) 
{
  if (!is.numeric(x)) {
    stop("argument x is not numeric")
  }
  if (!is.numeric(y)) {
    stop("argument y is not numeric")
  }
  if (is.matrix(x)) {
    Xmat <- x
  }
  else {
    if (is.vector(x)) {
      Xmat <- matrix(x, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (is.matrix(y)) {
    Ymat <- y
  }
  else {
    if (is.vector(y)) {
      Ymat <- matrix(y, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (nrow(Xmat) != nrow(Ymat)) 
    stop("argumentx x and y do not have the same row order")
  if (ncol(Xmat) != ncol(Ymat)) 
    stop("arguments x and y do not have the same column order")
  return(Xmat * Ymat)
}

calcPredInt <- function(x, sigma2, tau2, k, level = 0.95) {
  pi.lb <- x + qt((1-level)/2, df=(k-2))*sqrt(tau2+sigma2)
  pi.ub <- x + qt((1+level)/2, df=(k-2))*sqrt(tau2+sigma2)
  out <- list(lower = as.numeric(pi.lb), upper = as.numeric(pi.ub))
  return(out)
}




generateHyperparametersMA <- function(x, ...) {
  # Generate initial values from the relevant distributions
  model.pars <- list()
  model.pars[[1]] <- list(param = "mu.tobs", 
                          param.f = rnorm, 
                          param.args = list(n = 1, 
                                            mean = x$hp.mu.mean, 
                                            sd = sqrt(x$hp.mu.var)))
  
  if (x$hp.tau.dist == "dunif") {
    model.pars[[2]] <- list(param = "bsTau", 
                            param.f = runif, 
                            param.args = list(n = 1, 
                                              min = x$hp.tau.min, 
                                              max = x$hp.tau.max))
  } else if (x$hp.tau.dist == "dhalft") {
    model.pars[[2]] <- list(param = "bsTau", 
                            param.f = rstudentt, 
                            param.args = list(n = 1, 
                                              mean = x$hp.tau.mean, 
                                              sigma = x$hp.tau.sigma, 
                                              df = x$hp.tau.df,
                                              lower = x$hp.tau.min, 
                                              upper = x$hp.tau.max))
  } else {
    stop("Invalid distribution for 'hp.tau.dist'!")
  }
  model.pars
}


restore.c.var.hanley <- function(cstat, N.subjects, N.events, restore.method=4, g=NULL) {
  
  fHanley <- function(cstat, nstar, mstar, m, n) {
    ((cstat*(1-cstat)*(1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n))
  }
  
  n <- N.events #Number of events
  m <- N.subjects-N.events #Number of non-events
  
  if (restore.method==2) {
    mstar <- m-1
    nstar <- n-1
  } else if (restore.method==4) {
    mstar <- nstar <- N.subjects/2-1
  } else {
    stop ("Method not implemented yet!")
  }
  
  # Crude estimate
  cstat.var <- fHanley(cstat=cstat, m=m, n=n, mstar=mstar, nstar=nstar)
  if (is.null(g)) {
    return(cstat.var)
  }
  
  # Apply delta method if a transformation for c is defined
  ti <- rep(NA, length(cstat))
  for (i in 1:length(cstat)) {
    ci <- cstat[i]
    vi <- cstat.var[i]
    names(ci) <- names(vi) <- "cstat"
    ti[i] <- as.numeric((deltaMethod(object=ci, g=g, vcov.=vi))["SE"])**2
  }
  return(ti)
}


restore.c.var.se <- function(cstat, c.se, g=NULL) {
  if(is.null(g)) {
    return (c.se**2)
  }
  
  ti <- rep(NA, length(cstat))
  
  for (i in 1:length(cstat)) {
    ci <- cstat[i]
    vi <- c.se[i]**2
    names(ci) <- names(vi) <- "cstat"
    ti[i] <- as.numeric((deltaMethod(object=ci, g=g, vcov.=vi))["SE"])**2
  }
  return(ti)
}

restore.c.var.hanley2<- function(sd.LP, N.subjects, N.events, restore.method=4, g=NULL) {
  cstat <- calculate.cstat.sdPI(sd.LP, g=g)
  return(restore.c.var.hanley(cstat=cstat, N.subjects=N.subjects, N.events=N.events, restore.method=restore.method, g=g))
}

restore.c.var.se <- function(cstat, c.se, g=NULL) {
  if(is.null(g)) {
    return (c.se**2)
  }
  
  ti <- rep(NA, length(cstat))
  
  for (i in 1:length(cstat)) {
    ci <- cstat[i]
    vi <- c.se[i]**2
    names(ci) <- names(vi) <- "cstat"
    ti[i] <- as.numeric((deltaMethod(object=ci, g=g, vcov.=vi))["SE"])**2
  }
  return(ti)
}

restore.c.var.ci <- function(cil, ciu, level, g=NULL) {
  if(!is.null(g)) {
    lower <- eval(parse(text=g), list(cstat = cil))
    upper <- eval(parse(text=g), list(cstat = ciu))
  } else {
    lower <- cil
    upper <- ciu
  }
  if(missing(level)) level <- rep(0.95, length(cil))
  
  return(((upper - lower)/(2*qnorm((1-level)/2)))**2)
}

calculate.cstat.theta <- function(cstat, g=NULL) {
  if(is.null(g)) {
    return(cstat)
  }
  return(eval(parse(text=g), list(cstat = cstat)))
}

calculate.cstat.sdPI <- function (sdPI, g=NULL) {
  myfun <- function(x, sd.lp) {
    inv.logit(sqrt(2)*sd.lp*x)*dnorm(x, mean=0, sd=1)
  }
  cstat <- rep(NA, length(sdPI))
  for (i in 1:length(sdPI)) {
    cstat[i] <- ifelse(is.na(sdPI[i]), NA, 2*(integrate(myfun, lower=0, upper=+Inf, sd.lp=sdPI[i]))$value)
  }
  if(is.null(g)) {
    return(cstat)
  }
  return(eval(parse(text=g), list(cstat = cstat)))
}


# See params of geom_smooth for more details
# Smoothed calibration plot: use  formula = obsy ~ splines::bs(predy, 3)
plotCalibration <- function(predy, obsy, modelname="Model", 
                            formula = obsy ~ predy, 
                            method="glm", se = se, 
                            level=0.95, 
                            fam=binomial, ...) {
  #require(ggplot2)
  #require(ggExtra) ## Do the extra outside the package
  
  # Use formula instead
  
  # Add density plot under gg plot
  #predy<-rnorm(300)
  #obsy<-rt(300,df=10)
  #family <- gaussian
  
  xy <- data.frame(predy, obsy)
  
  scatter <- ggplot(xy, aes(x=predy, y=obsy)) +
    labs(x = "Predicted", y="Observed") + 
    scale_x_continuous(limits=c(min(predy),max(obsy))) + 
    scale_y_continuous(limits=c(min(predy),max(obsy))) #+ 
  
  
  # Add reference line for perfect calibration
  scatter <- scatter + geom_abline(aes(slope=1, intercept=0, linetype="Perfect calibration"), size=1)
  scatter <- scatter + geom_smooth(aes(linetype=modelname), method = method, 
                                   se = se, level=level, method.args = list(family = fam), ...)

  
  #scatter <- scatter + ggMarginal(data = xy, x = "predy", y = "obsy", margins = "x", type="histogram", size=4)
  #scatter <- scatter + labs(x = "Predicted", y="Observed") 
  scatter
}

