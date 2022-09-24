generateMCMCinits <- function(n.chains, model.pars, seed)
{
  inits <- list()
  for (i in seq_len(n.chains)) {
    inits.i <- list()
    for (j in seq_along(model.pars)) {
      parname <- model.pars[[j]]$param
      fprior <- model.pars[[j]]$param.f
      fargs <- model.pars[[j]]$param.args
      inits.i[[parname]] <- do.call(fprior, fargs)
    }
    
    # Do we need to add seed info?
    if (!is.null(seed)) {
      inits.i[[".RNG.name"]] <- "base::Wichmann-Hill"
      inits.i[[".RNG.seed"]] <- (seed + i) # Make sure each chain has a unique seed
    }
    
    inits[[i]] <- inits.i
  }
  
  return(inits)
}

.initiateDefaultPars <- function(pars, type = "") {
  pars.default <- list(level = 0.95,
                       hp.mu.mean = 0, 
                       hp.mu.var = 1000,
                       hp.tau.min = 0,
                       hp.tau.max = 100,
                       hp.tau.mean = 0,
                       hp.tau.sigma = 0.5,
                       hp.tau.dist = "dunif", 
                       hp.tau.df = 3, 
                       correction = 0.5,
                       seed = NULL)
  
  if (type == "valmeta") {
    pars.default$hp.tau.max = 2
    pars.default$method.restore.c.se = 4
    pars.default$model.cstat = "normal/logit"
    pars.default$model.oe = "normal/log"  #Alternative: "poisson/log" or "normal/identity"
  }
  
  if (!missing(pars)) {
    for (i in seq_along(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[element]]
    }
  }
  
  # Set random seed
  if (!is.null(pars.default$seed)) {
    message("Setting random seed to ", pars.default$seed)
    set.seed(pars.default$seed)
  }
  
  if (pars.default$level < 0 | pars.default$level > 1) {
    stop("Invalid value for 'level'!")
  } 
  
  return(pars.default)
}

# Random effects meta-analysis of regression coefficients
run_Bayesian_REMA <- function(call, measure, method, data, pars, n.chains, verbose, 
                              FUN_generate_bugs, # Function to generate BUGS code
                              ...) {
  
  # Identify number of studies
  numstudies = length(data$theta)
  
  # Perform a Bayesian meta-analysis
  model <- FUN_generate_bugs(pars = pars, ...)
  
  # Construct the hyperparameters
  model.pars <- generateHyperparametersMA(pars)
  
  # Generate initial values from the relevant distributions
  inits <- generateMCMCinits(n.chains = n.chains, 
                             model.pars = model.pars,
                             seed = pars$seed)
  
  # Generate data fame
  bugs_data <-  list(theta = data$theta,
                     theta.var = data$theta.se**2,
                     Nstudies = numstudies)
  
  jags.model <- runjags::run.jags(model = model$model.text, 
                                  monitor = c(model$model.pars, "PED"), 
                                  data = bugs_data, 
                                  confidence =  pars$level , # Which credibility intervals do we need?
                                  n.chains = n.chains,
                                  silent.jags = !verbose,
                                  inits = inits,
                                  ...)
  
  # Check if model converged
  psrf.ul <-  jags.model$psrf$psrf[,2]
  psrf.target <- jags.model$psrf$psrf.target
  
  if(sum(psrf.ul > psrf.target)>1) {
    warning(paste("Model did not properly converge! The upper bound of the convergence diagnostic (psrf) exceeds", 
                  psrf.target, "for the parameters", 
                  paste(rownames(jags.model$psrf$psrf)[which(psrf.ul > psrf.target)], " (psrf=", 
                        round(jags.model$psrf$psrf[which(psrf.ul > psrf.target),2],2), ")", collapse=", ", sep=""),
                  ". Consider re-running the analysis by increasing the optional arguments 'adapt', 'burnin' and/or 'sample'.", sep=""))
  }
  
  fit <- jags.model$summaries
  
  #Extract PED
  fit.dev <- runjags::extract(jags.model,"PED")
  txtLevel <- (pars$level*100)

  
  out <- list(call = call,
              model = pars$model.cstat,
              level = pars$level,
              numstudies = numstudies, 
              fit = jags.model, 
              PED = sum(fit.dev$deviance) + sum(fit.dev$penalty),
              est = fit[model$model.pars["mu_t"], "Mean"],
              se = fit[model$model.pars["mu_t"], "SD"],
              tau2 = fit[model$model.pars["tau2"], "Mean"],
              se.tau2 = fit[model$model.pars["tau2"], "SD"],
              ci.lb  = fit[model$model.pars["mu_t"], paste("Lower", txtLevel, sep = "")],
              ci.ub  = fit[model$model.pars["mu_t"], paste("Upper", txtLevel, sep = "")],
              pi.lb  = fit[model$model.pars["theta_new_t"],  paste("Lower", txtLevel, sep = "")],
              pi.ub  = fit[model$model.pars["theta_new_t"],  paste("Upper", txtLevel, sep = "")],
              data = data)
  
  class(out) <- "uvmeta"
  
  if (!missing(measure)) {
    out$measure <- measure
    
    if (measure == "cstat") {
      class(out) <- "valmeta"
    }
  }
  if (!missing(method)) {
    out$method <- method
  }
  
  return(out)
}

run_Bayesian_MA_oe <- function(call, measure, method, data, pars, n.chains, verbose, ...) {
  
  # Truncate hyper parameter variance
  pars$hp.mu.var = min(pars$hp.mu.var, 100)
  
  # Select studies where we have info on O, E and N
  i.select1 <- which(!is.na(data$O) & !is.na(data$E) & !is.na(data$N))
  
  # Select studies where we only have info on O and E
  i.select2 <- which(!is.na(data$O) & !is.na(data$E) & is.na(data$N))
  
  # Select studies where we have (estimated) information on log(OE) and its standard error
  i.select3 <- which(!is.na(data$theta) & !is.na(data$theta.se) & is.na(data$O) & is.na(data$E))
  
  mvmeta_dat <- list(O = data$O, E = data$E)
  
  if (length(i.select1) > 0) {
    mvmeta_dat$s1 <- i.select1
    mvmeta_dat$N <- data$N
  }
  if (length(i.select2) > 0)
    mvmeta_dat$s2 <- i.select2
  if (length(i.select3) > 0) {
    mvmeta_dat$s3 <- i.select3
    mvmeta_dat$logOE <- data$theta
    mvmeta_dat$logOE.se <- data$theta.se
  }
  
  # Generate model
  model <- generateBUGS.OE.discrete(N.type1 = length(i.select1), 
                                    N.type2 = length(i.select2),
                                    N.type3 = length(i.select3),
                                    pars = pars, ...)
  
  # Generate initial values from the relevant distributions
  model.pars <- generateHyperparametersMA(pars, ...)
  
  
  inits <- generateMCMCinits(n.chains = n.chains, 
                             model.pars = model.pars,
                             seed = pars$seed)
  
  jags.model <- runjags::run.jags(model = model, 
                                  monitor = c("mu.tobs", "mu.oe", "pred.oe", "bsTau", "prior_bsTau", "prior_mu", "PED"), 
                                  data = mvmeta_dat, 
                                  n.chains = n.chains,
                                  confidence = pars$level, # Which credibility intervals do we need?
                                  silent.jags = !verbose,
                                  inits = inits,
                                  ...)
  # cat(paste("\nPenalized expected deviance: ", round(x$PED,2), "\n"))
  
  # Check convergence
  psrf.ul <-  jags.model$psrf$psrf[,2]
  psrf.target <- jags.model$psrf$psrf.target
  
  if (sum(psrf.ul > psrf.target) > 0) {
    warning(paste("Model did not properly converge! The upper bound of the convergence diagnostic (psrf) exceeds", 
                  psrf.target, "for the parameters", 
                  paste(rownames(jags.model$psrf$psrf)[which(psrf.ul > psrf.target)], " (psrf=", 
                        round(jags.model$psrf$psrf[which(psrf.ul > psrf.target),2],2), ")", collapse = ", ", sep = ""),
                  ". Consider re-running the analysis by increasing the optional arguments 'adapt', 'burnin' and/or 'sample'."  ))
  }
  
  fit <- jags.model$summaries
  
  #Extract PED
  fit.dev <- runjags::extract(jags.model,"PED")
  txtLevel <- (pars$level*100)
  
  out <- list(call = call,
              measure = measure,
              method = method,
              model = pars$model.oe,
              level = pars$level,
              numstudies = length(c(i.select1, i.select2, i.select3)), 
              fit = jags.model, 
              PED = sum(fit.dev$deviance) + sum(fit.dev$penalty),
              est = fit["mu.oe", "Median"],
              ci.lb  = fit["mu.oe", paste("Lower", txtLevel, sep = "")],
              ci.ub  = fit["mu.oe", paste("Upper", txtLevel, sep = "")],
              pi.lb  = fit["pred.oe", paste("Lower", txtLevel, sep = "")],
              pi.ub  = fit["pred.oe", paste("Upper", txtLevel, sep = "")],
              data = data)
  class(out) <- "valmeta"
  return(out)
}

