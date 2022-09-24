.generateBugsCstat <- function(pars, 
                               ...) # standard deviation for student T prior
{
  
  hp.tau.prec <- 1/(pars$hp.tau.sigma**2)
  hp.mu.prec <- 1/pars$hp.mu.var
  
  out <- "model {\n " 
  out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
  out <- paste(out, "    theta[i] ~ dnorm(alpha[i], wsprec[i])\n")
  out <- paste(out, "    alpha[i] ~ dnorm(mu.tobs, bsprec)\n")
  out <- paste(out, "    wsprec[i] <- 1/(theta.var[i])\n")
  out <- paste(out, " }\n")
  out <- paste(out, " bsprec <- 1/(bsTau*bsTau)\n")
  out <- paste(out, " bsTauSq <- bsTau*bsTau\n")
  
  if (pars$hp.tau.dist == "dunif") {
    out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
    out <- paste(out, "  prior_bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "")
  } else if (pars$hp.tau.dist == "dhalft") {
    out <- paste(out, "  bsTau ~ dt(", pars$hp.tau.mean," ,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
    out <- paste(out, "  prior_bsTau ~ dt(", pars$hp.tau.mean," ,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "")
  } else {
    stop("Specified prior not implemented")
  }
  
  if (pars$model.cstat  == "normal/logit") {
    out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
    out <- paste(out, "  prior_mu ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
    out <- paste(out, "  mu.obs <- 1/(1+exp(-mu.tobs))\n", sep = "")
    out <- paste(out, "  pred.obs <- 1/(1+exp(-pred.tobs))\n", sep = "")
    out <- paste(out, "  pred.tobs ~ dnorm(mu.tobs, bsprec)\n", sep = "")
    
    model.pars <- c(mu = "mu.tobs", # Meta-analysis mean
                    mu_t = "mu.obs", # Transformed meta-analysis mean
                    tau2 = "bsTauSq", 
                    tau = "bsTau", 
                    prior_tau = "prior_bsTau", 
                    prior_mu = "prior_mu", 
                    theta_new = "pred.tobs", # New draw from the meta-analysis distribution
                    theta_new_t = "pred.obs" # New draw on the transformed scale
                    )
  } else if (pars$model.cstat == "normal/identity") {
    out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
    out <- paste(out, "  prior_mu ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
    out <- paste(out, "  pred.tobs ~ dnorm(mu.tobs, bsprec)\n", sep = "")
    
    model.pars <- c(mu = "mu.tobs", # Meta-analysis mean
                   mu_t = "mu.tobs", # Transformed meta-analysis mean
                   tau2 = "bsTauSq", 
                   tau = "bsTau", 
                   prior_tau = "prior_bsTau", 
                   prior_mu = "prior_mu", 
                   theta_new = "pred.tobs", # New draw from the meta-analysis distribution
                   theta_new_t = "pred.tobs" # New draw on the transformed scale
                  )
  } else {
    stop("Specified link function not implemented")
  }
  out <- paste(out, "}", sep = "")
  
  return(list(model.text = out, model.pars = model.pars))
}



generateBUGS.OE.discrete <- function(N.type1, N.type2, N.type3, pars, ...) {
  hp.tau.prec <- 1/(pars$hp.tau.sigma**2)
  hp.mu.prec <- 1/pars$hp.mu.var
  
  out <- "model {\n " 
  
  # Likelihood of studies providing O, E and N
  if (N.type1 > 0) {
    out <- paste(out, "for (j in 1:", N.type1, ")\n  {\n", sep="")
    out <- paste(out, "    O[s1[j]] ~ dbinom(pobs[j], N[s1[j]])\n")
    out <- paste(out, "    OE[j] <- exp(theta[s1[j]])\n")
    
    # Make sure that 'pobs' is always between 0 and 1
    out <- paste(out, "    pobs[j] <- min(OE[j], N[s1[j]]/(E[s1[j]]+1)) * E[s1[j]]/N[s1[j]] \n")
    out <- paste(out, " }\n")
  }
  
  # Likelihood of studies providing O and E
  if (N.type2 > 0) {
    out <- paste(out, "for (j in 1:", N.type2, ")\n  {\n", sep="")
    out <- paste(out, "    O[s2[j]] ~ dpois(lambda[j]) \n")
    out <- paste(out, "    lambda[j] <- exp(theta[s2[j]])*E[s2[j]]\n")
    out <- paste(out, " }\n")
  }
  
  # Likelihood for studies providing log OE ratio
  if (N.type3>0) {
    out <- paste(out, "for (j in 1:", N.type3, ")\n  {\n", sep="")
    out <- paste(out, "    logOE[s3[j]] ~ dnorm(theta[s3[j]], wsprec[j])\n")
    out <- paste(out, "    wsprec[j] <- 1/(logOE.se[s3[j]]*logOE.se[s3[j]])\n")
    out <- paste(out, " }\n")
  }
  
  # Between-study distribution of the logOR
  out <- paste(out, "for (j in 1:", (N.type1+N.type2+N.type3), ")\n  {\n", sep="")
  out <- paste(out, "    theta[j] ~ dnorm(mu.tobs, bsprec.logoe)\n")
  out <- paste(out, " }\n")
  out <- paste(out, " bsprec.logoe <- 1/(bsTau*bsTau)\n")
  if (pars$hp.tau.dist == "dunif") {
    out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
    out <- paste(out, "  prior_bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
  } else if (pars$hp.tau.dist == "dhalft") {
    out <- paste(out, "  bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
    out <- paste(out, "  prior_bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
  } else {
    stop("Specified prior not implemented")
  }
  out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
  out <- paste(out, "  prior_mu ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
  out <- paste(out, "  mu.oe <- exp(mu.tobs)\n", sep = "")
  out <- paste(out, "  pred.oe <- exp(pred.logoe)\n", sep = "")
  out <- paste(out, "  pred.logoe ~ dnorm(mu.tobs, bsprec.logoe)\n", sep = "")
  out <- paste(out, "}", sep = "")
  return(out)
}


generateBugsOE <- function(extrapolate=F,
                           pars,...) {
  hp.tau.prec <- 1/(pars$hp.tau.sigma**2)
  hp.mu.prec <- 1/pars$hp.mu.var
  
  out <- "model {\n " 
  
  if (pars$model.oe == "normal/log") {
    out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
    out <- paste(out, "    theta[i] ~ dnorm(alpha[i], wsprec[i])\n")
    out <- paste(out, "    alpha[i] ~ dnorm(mu.tobs, bsprec)\n")
    out <- paste(out, "    wsprec[i] <- 1/(theta.var[i])\n")
    out <- paste(out, " }\n")
    out <- paste(out, " bsprec <- 1/(bsTau*bsTau)\n")
    if (pars$hp.tau.dist == "dunif") {
      out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
      out <- paste(out, "  prior_bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "")
    } else if (pars$hp.tau.dist == "dhalft") {
      out <- paste(out, "  bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
      out <- paste(out, "  prior_bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
    } else {
      stop("Specified prior not implemented")
    }
    out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
    out <- paste(out, "  prior_mu ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
    out <- paste(out, "  mu.oe <- exp(mu.tobs)\n", sep="")
    out <- paste(out, "  pred.oe <- exp(pred.logoe)\n", sep="")
    out <- paste(out, "  pred.logoe ~ dnorm(mu.tobs, bsprec)\n", sep="")
    
  } else if (pars$model.oe=="poisson/log") {
    out <- paste(out, "for (i in 1:Nstudies)\n  {\n")
    out <- paste(out, "    obs[i] ~ dpois(mu[i])\n")
    out <- paste(out, "    mu[i] <- exc[i] * theta[i]\n")
    #out <- paste(out, "    exc[i] ~ dbinom (q[i], N[i])\n")
    #out <- paste(out, "    logit(q[i]) <- alphaQ[i]\n")
    #out <- paste(out, "    alphaQ[i] ~ dnorm(0.0,1.0E-6)\n")
    out <- paste(out, "    theta[i] <- exp(alpha[i])\n")
    out <- paste(out, "    alpha[i] ~ dnorm(mu.tobs, bsprec)\n")
    out <- paste(out, " }\n")
    out <- paste(out, " bsprec <- 1/(bsTau*bsTau)\n")
    if (pars$hp.tau.dist=="dunif") {
      out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
      out <- paste(out, "  prior_bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
    } else if (pars$hp.tau.dist=="dhalft") {
      out <- paste(out, "  bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "")
      out <- paste(out, "  prior_bsTau ~ dt(0,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
    } else {
      stop("Specified prior not implemented")
    }
    
    out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
    out <- paste(out, "  prior_mu ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
    out <- paste(out, "  mu.oe <- exp(mu.tobs)\n", sep = "")
    out <- paste(out, "  pred.oe <- exp(pred.logoe)\n", sep = "")
    out <- paste(out, "  pred.logoe ~ dnorm(mu.tobs, bsprec)\n", sep = "")
  }
  
  out <- paste(out, "}", sep = "")
  return(out)
}


.generateBugsREMA <- function(pars, 
                               ...) # standard deviation for student T prior
{
  hp.tau.prec <- 1/(pars$hp.tau.sigma**2)
  hp.mu.prec <- 1/pars$hp.mu.var
  
  out <- "model {\n " 
  out <- paste(out, "for (i in 1:Nstudies) {\n")
  out <- paste(out, "    wsprec[i] <- 1/theta.var[i]\n")
  out <- paste(out, "    theta[i] ~ dnorm(alpha[i], wsprec[i])\n")
  out <- paste(out, "    alpha[i] ~ dnorm(mu.tobs,prec)\n")
  out <- paste(out, " }\n\n")
  out <- paste(out, " #prior distributions\n")
  out <- paste(out, " tausq <- bsTau*bsTau\n")
  out <- paste(out, " prec <- 1/(tausq)\n")
  
  if (pars$hp.tau.dist == "dunif") {
    out <- paste(out, "  bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
    out <- paste(out, "  prior_bsTau ~ dunif(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "")
  } else if (pars$hp.tau.dist == "dhalft") {
    out <- paste(out, "  bsTau ~ dt(", pars$hp.tau.mean," ,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "") 
    out <- paste(out, "  prior_bsTau ~ dt(", pars$hp.tau.mean," ,", hp.tau.prec, ",", pars$hp.tau.df, ")T(", pars$hp.tau.min, ",", pars$hp.tau.max, ")\n", sep = "")
  } else {
    stop("Specified prior not implemented")
  }
  out <- paste(out, "  mu.tobs ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")
  out <- paste(out, "  prior_mu ~ dnorm(", pars$hp.mu.mean, ",", hp.mu.prec, ")\n", sep = "")

  out <- paste(out, " # predictive distribution\n")
  out <- paste(out, "  theta.new  ~ dnorm(mu.tobs, prec)\n", sep = "")

  out <- paste(out, "}", sep = "")
  
  ret.out <- list(model.text = out, 
                  model.pars = c(mu = "mu.tobs", # Meta-analysis mean
                                 mu_t = "mu.tobs", # Transformed meta-analysis mean
                                 tau2 = "tausq", 
                                 tau = "bsTau", 
                                 prior_tau = "prior_bsTau", 
                                 prior_mu = "prior_mu", 
                                 theta_new = "theta.new", # New draw from the meta-analysis distribution
                                 theta_new_t = "theta.new" # New draw on the transformed scale
                                 )
                  )
  
  return(ret.out)
}