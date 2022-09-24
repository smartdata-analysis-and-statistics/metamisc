#' Impute missing values by their conditional mean
#'
#' This function imputes missing values by their conditional mean
#' 
#' @param x A vector with observations, some of which may be missing (indicated by NA)
#' @param mu A vector with the population means for 'x'. No missing values are allowed here.
#' @param Sigma A matrix describing the population covariance of 'x'
#' 
#' @return A vector where missing values for 'x' have been replaced by their conditional mean
#' 
#' @examples 
#' # Define the population means
#' mu <- c(0, 1, 2)
#' 
#' # Define the covariance of the population
#' Sigma <- diag(1,3)
#' Sigma[1,2] <- Sigma[2,1] <- 0.3 
#' Sigma[2,3] <- Sigma[3,2] <- 0.1
#' Sigma[1,3] <- Sigma[3,1] <- -0.2
#' 
#' # Generate a 'random' sample from the population that is partially observed
#' x <- c(NA, 2, 4)
#' 
#' # Impute the missing values
#' impute_conditional_mean (x=x, mu=mu, Sigma=Sigma)
#' 
#' @keywords imputation
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
impute_conditional_mean <- function(x, mu, Sigma) {
  
  
  if (!is.array(Sigma) & !is.matrix(Sigma)) {
    stop ("Sigma should be a matrix!")
  }
  
  # Verify that x and mu are a vector
  if(!(is.vector(x) & is.vector(mu))) {
    stop ("'x' and 'mu' should be a vector!")
  }
  
  # Check the dimension of all input variables
  if (length(unique(c(length(x), length(mu), nrow(sigma), ncol(sigma))))>1) {
    stop ("Incompatible size of input variables!")
  }
  
  # If x, mu and Sigma are named, make sure that they are in the same order
  if(!is.null(names(x)) & !is.null(names(mu))) {
    stopifnot(all.equal(names(x), names(mu)))
  }
  
  if (!is.null(names(mu)) & !is.null(colnames(Sigma))) {
    stopifnot(all.equal(names(mu), colnames(Sigma)))
    stopifnot(all.equal(names(mu), rownames(Sigma)))
  }
  
  if (!is.null(names(x)) & !is.null(colnames(Sigma))) {
    stopifnot(all.equal(names(x), colnames(Sigma)))
    stopifnot(all.equal(names(x), rownames(Sigma)))
  }
  
  # Return 'x' if no imputation is required
  if(sum(is.na(x))==0) {
    return (x)
  }
  
  # If all values for 'x' are missing, simply return the means
  if (sum(is.na(x))==length(x)) {
    return (mu)
  }
  
  # Check if Sigma is symmetric
  if (!isSymmetric(Sigma)) {
    stop ("The covariance matrix is not symmetric!")
  }
  
  # Check if the covariance matrix is valid
  if (any(eigen(Sigma, only.values = TRUE)$values <= 0 )) {
    stop ("The covariance matrix is not positive semidefinite!")
  }
  
  # Identify which values of x are missing and should be imputed
  dependent_ind <- which(is.na(x))
  given_ind <- which(!is.na(x))
  
  C <- Sigma[dependent_ind, given_ind, drop = FALSE]
  D <- Sigma[given_ind, given_ind]
  CDinv <- C %*% solve(D)
  cMu <- c(mu[dependent_ind] + CDinv %*% (x[given_ind] - mu[given_ind]))
  
  # Impute
  x_imputed <- x
  x_imputed[dependent_ind] <- cMu
  x_imputed
}
 