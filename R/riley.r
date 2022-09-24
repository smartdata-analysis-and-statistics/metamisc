#TODO: check prediction intervals
#TODO: Add cholesky decomposition in loglikelihood
#TODO: Alter design matrix in LogLik to set entries with missing data to zero (in which case #df needs to be altered)
#TODO: Calculate SE of rho using Delta method


#' Fit the alternative model for bivariate random-effects meta-analysis
#' 
#' This function fits the alternative model for bivariate random-effects meta-analysis when the within-study correlations 
#' are unknown. This bivariate model was proposed by Riley et al. (2008) and is similar to the general bivariate 
#' random-effects model (van Houwelingen et al. 2002), but includes an overall correlation parameter rather than 
#' separating the (usually unknown) within- and between-study correlation. As a consequence, the alternative model 
#' is not fully hierarchical, and estimates of additional variation beyond sampling error (\code{psi}) are not 
#' directly equivalent to the between-study variation (\code{tau}) from the general model. This model is particularly 
#' useful when there is large within-study variability, few primary studies are available or the general model 
#' estimates the between-study correlation as 1 or -1. 
#' 
#' @param  X data frame containing integer variables \code{Y1}, \code{vars1}, \code{Y2} and \code{vars2},
#' where the columns \code{Y1} and \code{Y2} represent the effect sizes of outcome 1 and, respectively, outcome 2. The columns 
#' \code{vars1} and \code{vars2} represent the error variances of \code{Y1} and, respectively, \code{Y2}. Alternatively, 
#' when considering a meta-analysis of diagnostic test accuracy data, the columns \code{TP}, \code{FN}, \code{FP} and 
#' \code{TN} may be specified. Corresponding values then represent the number of true positives, the number of false negatives,
#' the number of false positives and, respectively, the number of true negatives.
#' @param slab Optional vector specifying the label for each study
#' @param optimization The optimization method that should be used for minimizing the negative (restricted) 
#' log-likelihood function. The default method is an implementation of that of Nelder and Mead (1965), 
#' that uses only function values and is robust but relatively slow. Other methods are described in \link[stats]{optim}.
#' @param control A list of control parameters to pass to \link[stats]{optim}.
#' @param pars List with additional arguments. The width of confidence, credibility and prediction intervals is 
#' defined by \code{level} (defaults to 0.95). 
#' @param \dots Arguments to be passed on to other functions. See "Details" for more information.
#' 
#' @details Parameters are estimated by iteratively maximizing the restriced log-likelihood using the Newton-Raphson procedure. 
#' The results from a univariate random-effects meta-analysis with a method-of-moments estimator are used as starting 
#' values for \code{beta1}, \code{beta2}, \code{psi1} and \code{psi2} in the \code{optim} command. Standard errors for all parameters are obtained 
#' from the inverse Hessian matrix.
#' 
#' \subsection{Meta-analysis of effect sizes}{
#' The following parameters are estimated by iteratively maximizing the restriced log-likelihood using the Newton-Raphson 
#' procedure: pooled effect size for outcome 1 (\code{beta1}), pooled effect size for outcome 2 (\code{beta2}), 
#' additional variation of \code{beta1} beyond sampling error (\code{psi1}), additional variation of \code{beta2} 
#' beyond sampling error (\code{psi2}) and the correlation \code{rho} between \code{psi1} and \code{psi2}. 
#' 
#' }
#' 
#' \subsection{Meta-analysis of diagnostic test accuracy}{
#' Although the model can also be used for diagnostic test accuracy data when substantial within-study correlations 
#' are expected, assuming zero within-study correlations (i.e. applying Reitsma's approach) is usually justified 
#' (Reitsma et al. 2005, Daniels and Hughes 1997, Korn et al. 2005, Thompson et al. 2005, Van Houwelingen et al. 2002).
#' 
#' A logit transformation is applied to the sensitivities ans false positive rates of each study, in order to meet the normality 
#' assumptions. When zero cell counts occur, continuity corrections may be required. The correction value can be specified using
#' \code{correction} (defaults to 0.5). Further, when the argument \code{correction.control} is set to \code{"all"} 
#' (the default) the continuity correction is added to the whole data if only one cell in one study is zero. 
#' If \code{correction.control="single"} the correction is only applied to rows of the data which have a zero.
#' 
#' The following parameters are estimated: logit of sensitivity (\code{beta1}), logit of false positive rate (\code{beta2}), 
#' additional variation of \code{beta1} beyond sampling error (\code{psi1}), additional variation of \code{beta2} beyond 
#' sampling error (\code{psi2}) and the correlation (\code{rho}) between \code{psi1} and \code{psi2}. 
#' }
#' 
#' @note The overall correlation parameter \code{rho} is transformed during estimation to ensure that corresponding values
#' remain between -1 and 1. The transformed correlation \code{rhoT} is given as \code{logit((rho+1)/2)}. During optimization,
#' the starting value for \code{rhoT} is set to 0.  The standard error of \code{rho} is derived from \code{rhoT} using 
#' the Delta method. Similarly, the Delta methods is used to derive the standard error of the sensitivity and false 
#' positive rate from \code{beta1} and, respectively, \code{beta2}.
#' 
#' Algorithms for dealing with missing data are currently not implemented, but Bayesian approaches will become 
#' available in later versions. 
#' 
#' @references
#' \itemize{
#' \item Korn EL, Albert PS, McShane LM. Assessing surrogates as trial endpoints using mixed models. 
#' \emph{Statistics in Medicine} 2005; \bold{24}: 163--182.
#' \item Nelder JA, Mead R. A simplex algorithm for function minimization. \emph{Computer Journal} (1965); \bold{7}: 308--313.
#' \item Reitsma J, Glas A, Rutjes A, Scholten R, Bossuyt P, Zwinderman A. Bivariate analysis of sensitivity and 
#' specificity produces informative summary measures in diagnostic reviews. \emph{Journal of Clinical Epidemiology} 2005; 
#' \bold{58}: 982--990.
#' \item Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when 
#' the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
#' \item Thompson JR, Minelli C, Abrams KR, Tobin MD, Riley RD. Meta-analysis of genetic studies using mendelian 
#' randomization--a multivariate approach. \emph{Statistics in Medicine} 2005; \bold{24}: 2241--2254.
#' \item van Houwelingen HC, Arends LR, Stijnen T. Advanced methods in meta-analysis: multivariate approach and 
#' meta-regression. \emph{Statistics in Medicine} 2002; \bold{21}: 589--624.
#' }
#' 
#' @examples 
#' data(Scheidler)
#' data(Daniels)
#' data(Kertai)
#' 
#' # Meta-analysis of potential surrogate markers data
#' # The results obtained by Riley (2008) were as follows:
#' # beta1 = -0.042 (SE = 0.063),
#' # beta2 = 14.072 (SE = 4.871)
#' # rho   = -0.759
#' \dontrun{
#' fit1 <- riley(Daniels) #maxit reached, try again with more iterations
#' }
#' fit1 <- riley(Daniels, control=list(maxit=10000))
#' summary(fit1)
#' 
#' # Meta-analysis of prognostic test studies
#' fit2 <- riley(Kertai)
#' fit2
#' 
#' # Meta-analysis of computed tomography data 
#' ds <- Scheidler[which(Scheidler$modality==1),]
#' fit3 <- riley(ds)
#' fit3
#' 
#' @return An object of the class \code{riley} for which many standard methods are available. A warning message is 
#' casted when the Hessian matrix contains negative eigenvalues, which implies that the identified solution is a 
#' saddle point and thus not optimal.
#' 
#' 
#' @keywords regression multivariate bivariate riley meta-analysis
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#'
#' @export
riley <- function(X, slab, optimization = "Nelder-Mead", control = list(), pars, ...) UseMethod("riley")

#' @export
riley.default <- function(X, slab, optimization = "Nelder-Mead", control = list(), pars, ...)
{
  out <- NA
  
  pars.default <- list(level = 0.95)
  
  # Load default parameters
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[element]]
    }
  }
  
  # Check parameter values
  if (pars.default$level < 0 | pars.default$level > 1) {
    stop ("Invalid value for significance level!")
  } 
  
  if(missing(X)) stop("X is missing!")
  
  k <- dim(X)[1] # Number of studies
  
  # Assess which type of meta-analysis is needed
  if (sum(c("Y1","vars1","Y2","vars2","TP","FN","FP","TN") %in% colnames(X))==8) {
    stop(paste("Too many variables specified in X! Please choose whether to perform a meta-analysis of effect sizes",
               "or a meta-analysis of diagnostic test accuracy data!"))
  }
  if (sum(c("Y1","vars1","Y2","vars2") %in% colnames(X))==4) {
    out <- rileyES(X, optimization = optimization, control=control, pars=pars.default, ...)
  } else if (sum(c("Y1","vars1","Y2","vars2") %in% colnames(X))>0) {
    stop ("Missing column(s) in X!")
  } else if (sum(c("TP","FN","FP","TN") %in% colnames(X))==4) {
    out <- rileyDA(X, optimization = optimization, control=control, pars=pars.default, ...)
  } else if (sum(c("TP","FN","FP","TN") %in% colnames(X))>0) {
    stop ("Missing column(s) in X!")
  } else {
    stop ("Provided data not supported, please verify column names in X!")
  }
  
  #######################################################################################
  # Assign study labels
  #######################################################################################
  if(missing(slab)) {
    out$slab <- paste("Study", seq(1, k))
  } else {
    out$slab <- make.unique(as.character(slab))
  }
  
  out$call <- match.call()
  
	class(out) <- "riley"
	out
}

# X is the design matrix
negfullloglikRiley <- function(parsll, X, Y,vars)
{

  Beta = rbind(parsll[1], parsll[2]) #Beta vector
  psisq1 = parsll[3]**2 #ensure variance is positive
  psisq2 = parsll[4]**2 #ensure variance is positive
  rho = inv.logit(parsll[5])*2-1 #ensure correlation is in [-1,1], and values in that interval move symmetric from -1 to 0 and from 1 to 0
  k = 2 #2 endpoints
  n = dim(Y)[1]/2
  
  #Create Phi matrix
  Phi = array(0,dim=c((n*k),(n*k)))
  for (i in 1:n) {
    Phi[((i-1)*2+1),((i-1)*2+1)] = vars[i,1]+psisq1
    Phi[((i-1)*2+2),((i-1)*2+2)] = vars[i,2]+psisq2
    Phi[((i-1)*2+1),((i-1)*2+2)] = rho*sqrt((vars[i,1]+psisq1)*(vars[i,2]+psisq2))
    Phi[((i-1)*2+2),((i-1)*2+1)] = rho*sqrt((vars[i,1]+psisq1)*(vars[i,2]+psisq2))
  }
  
  #Minimize the negative of the restricted log-lkh
  0.5*((n-k)*log(2*pi)-log(det(t(X)%*%X))+log(det(Phi))+log(det(t(X)%*%solve(Phi)%*%X))+(t(Y-X%*%Beta)%*%solve(Phi)%*%(Y-X%*%Beta)))
}

# effect sizes data meta-analysis
rileyES <- function(X = NULL, Y1, Y2, vars1, vars2, optimization = "Nelder-Mead", control = list(), pars=list(level=0.95), ...)
{
	if(!is.null(X)){
		X <- as.data.frame(X)
		origdata <- X
		Y1 <- X$Y1
		Y2 <- X$Y2
		vars1 <- X$vars1
		vars2 <- X$vars2
	} else {
		origdata <- cbind(Y1,vars1,Y2,vars2)
		colnames(origdata) <- c("Y1","vars1","Y2","vars2")
	}
  
	numstudies <- length(Y1)
	nobs <- length(which(!is.na(Y1)))+length(which(!is.na(Y2)))
	
	if(nobs != numstudies*2){warning("There are missing observations in the data!")}
	
	df <- 5 #There are 5 parameters to estimate
	if(numstudies*2-df < 0){warning("There are very few primary studies!")}
	
	# Set up the design matrix
	Xdesign <- array(0,dim=c(numstudies*2,2))
	Xdesign[seq(from=1, to=(numstudies*2), by=2)[which(!is.na(Y1))],1] <- 1 #Indicate which variables are observed for Y1
	Xdesign[seq(from=2, to=(numstudies*2), by=2)[which(!is.na(Y2))],2] <- 1
	
	
	vars <- cbind(vars1, vars2)
	Y <- array(NA,dim=c((length(Y1)*2),1))
	Y[seq(from=1, to=(numstudies*2), by=2)] <- Y1
	Y[seq(from=2, to=(numstudies*2), by=2)] <- Y2
	
	#Calculate starting values for optim
	pars.start = c(0,0,0,0,0)
	if (numstudies >= 2) {
		sumlY1 <- metafor::rma(yi=Y1, vi=vars1, method="DL")
		sumlY2 <- metafor::rma(yi=Y2, vi=vars2, method="DL")
		pars.start = c(sumlY1$beta[1], sumlY2$beta[1], sqrt(vcov(sumlY1)[1,1]), sqrt(vcov(sumlY2)[1,1]), 0)
	}
	
	
	
	fit = optim(pars.start, negfullloglikRiley, X=Xdesign, Y=Y, vars=vars, method=optimization, hessian=T, control=control)
	
	if(fit$convergence != 0) { 
		if(fit$convergence == 1) warning ("Iteration limit had been reached.")
		else if (fit$convergence == 10) warning("Degeneracy of the Nelder-Mead simplex.")
		else if (fit$convergence == 51 | fit$convergence == 52) warning(fit$message)
    else warning("Unspecified convergence error in optim.")
	}
	
	beta1 = fit$par[1]
	beta2 = fit$par[2]
	psi1 = abs(fit$par[3])
	psi2 = abs(fit$par[4])
	rhoT = fit$par[5]
	coefficients = c(beta1,beta2,psi1,psi2,rhoT)
	names(coefficients) = c("beta1","beta2","psi1","psi2","rhoT")
	
	hessian = fit$hessian
	colnames(hessian) = c("beta1","beta2","psi1","psi2","rhoT")
	rownames(hessian) = c("beta1","beta2","psi1","psi2","rhoT")
	
	if (length(which(eigen(fit$hessian,symmetric=TRUE)$values<0))>0) warning("The Hessian contains negative eigenvalues!")
	
	iterations <- fit$iterations
	logLik <- -fit$value
	
	output <- list(coefficients = coefficients, hessian = hessian, df = df, numstudies = numstudies, nobs = nobs, 
	               df.residual = (sum(Xdesign)-5), #Remainig degrees of freedom
	               logLik = logLik, iterations = (iterations+1), data = origdata, type="effect.size", level=pars$level)  
	return(output)
}

# Diagnostic test accuracy data meta-analysis
rileyDA <-
  function(X = NULL, TP, FN, FP, TN, correction = 0.5, 
           correction.control = "all", optimization = "Nelder-Mead", control = list(), pars=list(level=0.95), ...)
  {
      if(!is.null(X)){
        X <- as.data.frame(X)
        origdata <- newdata <- X
      } else {
        origdata <- newdata <- as.data.frame(cbind(TP,FN,FP,TN))
        colnames(origdata) <- c("TP","FN","FP","TN")
      }
      
	  ## The following corrections are copied from the "mada" package to facilitate comparison of results
      ## apply continuity correction to _all_ studies if one contains zero
      if(correction.control == "all"){if(any(origdata == 0)){newdata$TP <- origdata$TP + correction;
                                                             newdata$FN <- origdata$FN + correction;
                                                             newdata$FP <- origdata$FP + correction;
                                                             newdata$TN <- origdata$TN + correction}}
      if(correction.control == "single"){
        correction = ((((origdata$TP == 0)|(origdata$FN == 0))|(origdata$FP == 0))| (origdata$TN == 0))*correction
        newdata$TP <- correction + origdata$TP
        newdata$FN <- correction + origdata$FN
        newdata$FP <- correction + origdata$FP
        newdata$TN <- correction + origdata$TN
      }
      
      
      #Calculate sensitivities and specificities (original scale)
      number.of.pos <- newdata$TP + newdata$FN
      number.of.neg <- newdata$FP + newdata$TN
      sens <-newdata$TP/number.of.pos
      fpr <- newdata$FP/number.of.neg
      var.sens = sens*(1-sens)/number.of.pos
      var.fpr = fpr*(1-fpr)/number.of.neg
      
      logit.sens <- logit(sens)
      logit.fpr <- logit(fpr)
      var.logit.sens <- 1/(sens*(1-sens)*number.of.pos)
      var.logit.fpr <- 1/(fpr*(1-fpr)*number.of.neg)
      
	    #Apply ordinary bivariate meta-analysis on transformed data
      output = rileyES(X=NULL, Y1=logit.sens,Y2=logit.fpr,vars1=var.logit.sens,vars2=var.logit.fpr,optimization = optimization, control = control, ...)
      output$type = "test.accuracy"
      output$data = cbind(newdata, output$data)
      output$correction = correction 
      output$correction.control = correction.control
      output$level <- pars$level
      
      return(output)
  }

#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method print riley
#' @export
print.riley <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  # calculate rho and its standard error
  rho <- inv.logit(x$coefficients["rhoT"])*2-1
  rho.se <- estSErho(rhoT=x$coefficients["rhoT"], var.rhoT=sqrt(vcov(x)["rhoT","rhoT"]))
  
  cat("Call:\n")
  print(x$call)
  coefs <- x$coefficients
  coefs["rhoT"] <- inv.logit(x$coefficients["rhoT"])*2-1
  names(coefs)[which(names(coefs)=="rhoT")] <- "rho"
  cat("\nCoefficients\n")
  print(coefs)
  cat("\nDegrees of Freedom: ", x$df.residual, " Residual", sep="")
  
  if (length(which(eigen(x$hessian,symmetric=TRUE)$values<0))>0) cat("\nWarning: the Hessian matrix contains negative eigenvalues, parameter estimates are thus not optimally fitted!\n")
}


#' Parameter summaries
#' Provides the summary estimates of the alternative model for bivariate random-effects meta-analysis by Riley et al. 
#' (2008) with their corresponding standard errors (derived from the inverse Hessian). For confidence intervals,
#' asymptotic normality is assumed.
#' 
#' @param object A \code{riley} object
#' @param \dots Arguments to be passed on to other functions (currently ignored)
#' 
#' @details For meta-analysis of diagnostic test accuracy data, \code{beta1} equals the logit sensitivity (Sens) and 
#' \code{beta2} equals the logit false positive rate (FPR).
#' 
#' @note For the overall correlation (\code{rho}) confidence intervals are derived using the transformation 
#' \code{logit((rho+1)/2)}. Similarly, the logit transformation is used to derive confidence intervals for the summary
#' sensitivity and false positive rate.
#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method summary riley
#' 
#' @references 
#' Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the 
#' within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
#' 
#' @return array with confidence intervals for the estimated model parameters. For diagnostic test accuracy data, 
#' the resulting summary sensitivity and false positive rate are included.
#' 
#' @export
summary.riley <- function(object,  ...)
{
  # calculate rho and its standard error
  rho <- inv.logit(object$coefficients["rhoT"])*2-1
  rho.se <- estSErho(rhoT=object$coefficients["rhoT"], var.rhoT=sqrt(vcov(object)["rhoT","rhoT"]))
  
  results <- cbind(object$coefficients, sqrt(diag(vcov(object))), confint(object))
  colnames(results)[1:2] <- c("Estimate", "SE")
  
  # Change rhoT to rho
  results["rhoT", ] <- c(rho, rho.se, inv.logit(results["rhoT", 3:4])*2-1)
  rownames(results)[5] <- "rho"
  
  if (object$type=="test.accuracy") {
    se.sens <- estSEexpit(beta=object$coefficients["beta1"], var.beta=diag(vcov(object))["beta1"])
    se.fpr  <- estSEexpit(beta=object$coefficients["beta2"], var.beta=diag(vcov(object))["beta2"])
    results <- rbind(results, cbind(inv.logit(object$coefficients["beta1"]), se.sens, inv.logit(results["beta1", 3]), inv.logit(results["beta1", 4])))
    results <- rbind(results, cbind(inv.logit(object$coefficients["beta2"]), se.fpr, inv.logit(results["beta2", 3]), inv.logit(results["beta2", 4])))
    rownames(results)[6] <-  "Sens"
    rownames(results)[7] <-  "FPR"
  } 
 
  
  res <- list(call=object$call, confints = results)
  class(res) <- "summary.riley"
  res
}

estSEexpit <- function(beta, var.beta) {
  # Use Delta method to calculate SE(expit(beta))
  sqrt((((exp(beta)/(exp(beta)+1)**2))**2) * var.beta)
}

estSErho <- function(rhoT, var.rhoT) {
  # Use Delta method to calculate SE(rho)
  # rhoT = logit((rho+1)/2)
  # rho  = expit(rhoT)*2-1
  # var(rho) = (deriv(expit(rhoT)*2-1))**2 * var(rhoT)
  # 
  sqrt((((2*exp(rhoT)/(exp(rhoT)+1)**2))**2) * var.rhoT)
}

#' Prediction Interval
#' 
#' Calculates a prediction interval for the summary parameters of Riley's alternative model for bivariate random-effects 
#' meta-analysis. This interval predicts in what range future observations will fall given what has already been observed.
#' 
#' @param object  A \code{riley} object.
#' @param \dots Additional arguments (currently ignored)
#' 
#' @details Prediction intervals are based on Student's t-distribution with (numstudies - 5) degrees of freedom. The width of 
#' the interval is specified by the significance level chosen during meta-analysis.
#' 
#' @return Data frame containing prediction intervals with the summary estimates \code{beta1} and \code{beta2} 
#' (for effect size data), or with the mean sensitivity and false positive rate (for diagnostic test accuracy data).
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @method predict riley
#' @export
predict.riley <- function(object,  ...)
{
  alpha <- (1-object$level)/2
  
  predint		<- array(NA,dim=c(4,2))
  colnames(predint) <- c( paste((alpha*100),"%"),paste(((1-alpha)*100),"%"))
  df 			<- object$df
  numstudies 		<- object$numstudies
  
  rownames(predint) = c("beta1","beta2", "Sens","FPR")
  predint[1,] = c((qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"]))
  predint[2,] = c((qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"]))
  
  
  if (object$type=="test.accuracy") {
    if ((numstudies - df) > 0) {
      predint[3,] = inv.logit(c((qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"])))
      predint[4,] = inv.logit(c((qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"])))
    } else {
      predint[3,] = c(inv.logit(coefficients(object)["beta1"]),0,1)
      predint[4,] = c(inv.logit(coefficients(object)["beta2"]),0,1)
    }
  } else {
    predint <- predint[1:2,]
  }
  predint
}


print.summary.riley <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$confints)
}

#' @author Thomas Debray <thomas.debray@gmail.com>
#' @method vcov riley
#' @export
vcov.riley <- function(object, ...) {
  if (length(which(eigen(object$hessian,symmetric=TRUE)$values<0))>0) warning("The Hessian contains negative eigenvalues!")
 
  # It is known that 'optim' has problems.  Perhaps the simplest thing to do is to call 'optim' with each of 
  # the 'methods' in sequence, using the 'optim' found by each 'method' as the starting value for the next.  
  # When I do this, I often skip 'SANN', because it typically takes so much more time than the other methods.  
  # However, if there might be multiple local minima, then SANN may be the best way to find a global minimum, 
  # though you may want to call 'optim' again with another method, starting from optimal solution returned by 'SANN'. 

  Sigma = solve(object$hessian)
  Sigma
}

#' Print the log-likelihood
#' 
#' This function provides the (restricted) log-likelihood of a fitted model.
#' 
#' @param  object A \code{riley} object, representing a fitted alternative model for bivariate random-effects 
#' meta-analysis when the within-study correlations are unknown.
#' @param \dots Additional arguments to be passed on to other functions, currently ignored.
#' @return Returns an object of class \code{logLik}. This is the (restricted) log-likelihood of the model represented 
#' by \code{object} evaluated at the estimated coefficients. It contains at least one attribute, 
#' "\code{df}" (degrees of freedom), giving the number of (estimated) parameters in the model.
#' 
#' @references Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when 
#' the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
#' 
#' @examples 
#' data(Daniels)
#' fit <- riley(Daniels,control=list(maxit=10000))
#' logLik(fit)
#' 
#' @keywords likelihood riley bivariate meta-analysis
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @method logLik riley
#' @export
logLik.riley <- function(object, ...) {
	val 				      <- object$logLik
	attr(val, "nobs") <- object$nobs
	attr(val, "df") 	<- object$df
	class(val) 			  <- "logLik"
	return(val)
}

#' Plot the summary of the bivariate model from Riley et al. (2008).
#' 
#' Generates a forest plot for each outcome of the bivariate meta-analysis.
#' 
#' @param x An object of class \code{riley}
#' @param title Title of the forest plot
#' @param sort By default, studies are ordered by ascending effect size (\code{sort="asc"}). For study ordering by descending
#' effect size, choose \code{sort="desc"}. For any other value, study ordering is ignored.
#' @param xlim The \code{x} limits \code{(x1, x2)} of the forest plot
#' @param refline Optional numeric specifying a reference line
#' @param \dots Additional parameters for generating forest plots
#' 
#' @references Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis 
#' when the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
#' 
#' @examples 
#' data(Scheidler)
#' 
#' #Perform the analysis
#' fit <- riley(Scheidler[which(Scheidler$modality==1),])
#' plot(fit)
#' 
#' require(ggplot2)
#' plot(fit, sort="none", theme=theme_gray())
#' 
#' @keywords forest
#' @author Thomas Debray <thomas.debray@gmail.com>
#' @import grid
#' 
#' @method plot riley
#' @export
plot.riley <- function(x, title, sort="asc", xlim, refline, ...)
{
	alpha <- (1-x$level)/2
	
	if (x$type=="test.accuracy") {
	  if (missing(title)) {
	    title <- c("Sensitivity", "False Positive Rate")
	  }
	  if (missing(xlim)) {
	    xlim <- c(0, 1)
	  }

	  yi1    <- inv.logit(x$data[,"Y1"])
	  yi1.ci <- inv.logit(cbind(x$data[,"Y1"]+qnorm(alpha)*sqrt(x$data[,"vars1"]), x$data[,"Y1"]+qnorm(1-alpha)*sqrt(x$data[,"vars1"])))
	  ma1    <- inv.logit(coef(x)["beta1"])
	  ma1.ci <- inv.logit(confint(x)["beta1",])
	  ma1.pi <- inv.logit(predict(x)["beta1",])
	  
	  yi2    <- inv.logit(x$data[,"Y2"])
	  yi2.ci <- inv.logit(cbind(x$data[,"Y2"]+qnorm(alpha)*sqrt(x$data[,"vars2"]), x$data[,"Y2"]+qnorm(1-alpha)*sqrt(x$data[,"vars2"])))
	  ma2    <- inv.logit(coef(x)["beta2"])
	  ma2.ci <- inv.logit(confint(x)["beta2",])
	  ma2.pi <- inv.logit(predict(x)["beta2",])
	} else {
	  if (missing(title)) {
	    title <- c("Outcome 1", "Outcome 2")
	  }
	  if (missing(refline)) {
	    refline <- 0
	  }
	  
	  yi1    <- x$data[,"Y1"]
	  yi1.ci <- cbind(x$data[,"Y1"]+qnorm(alpha)*sqrt(x$data[,"vars1"]), x$data[,"Y1"]+qnorm(1-alpha)*sqrt(x$data[,"vars1"]))
	  ma1    <- coef(x)["beta1"]
	  ma1.ci <- confint(x)["beta1",]
	  ma1.pi <- predict(x)["beta1",]
	  
	  yi2    <- x$data[,"Y2"]
	  yi2.ci <- cbind(x$data[,"Y2"]+qnorm(alpha)*sqrt(x$data[,"vars2"]), x$data[,"Y2"]+qnorm(1-alpha)*sqrt(x$data[,"vars2"]))
	  ma2    <- coef(x)["beta2"]
	  ma2.ci <- confint(x)["beta2",]
	  ma2.pi <- predict(x)["beta2",]
	  
	}
	
	# Generate 2 forest plots
	p1 <- forest(theta=yi1, 
	             theta.ci.lb=yi1.ci[,1], theta.ci.ub=yi1.ci[,2], 
	             theta.slab=x$slab, 
	             theta.summary=ma1, 
	             theta.summary.ci.lb=ma1.ci[1], theta.summary.ci.ub=ma1.ci[2], 
	             theta.summary.pi.lb=ma1.pi[1], theta.summary.pi.ub=ma1.pi[2], 
	       title = title[1], xlim=xlim, refline=refline, sort=sort, ...)
	p2 <- forest(theta=yi2, 
	             theta.ci.lb=yi2.ci[,1], theta.ci.ub=yi2.ci[,2], 
	             theta.slab=x$slab, 
	       theta.summary=ma2, 
	       theta.summary.ci.lb=ma2.ci[1], theta.summary.ci.ub=ma2.ci[2], 
	       theta.summary.pi.lb=ma2.pi[1], theta.summary.pi.ub=ma2.pi[2], 
	       title = title[2], xlim=xlim, refline=refline, sort=sort, ...)
	
	multiplot(p1, p2, cols=2)
}
