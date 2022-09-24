\name{vcov.riley}
\alias{vcov.riley}

\title{Calculate Variance-Covariance Matrix for a Fitted Riley Model Object }

\description{ Returns the variance-covariance matrix of the main parameters of a fitted model object.}
\usage{\method{vcov}{riley}(object, \ldots)}

\arguments{
  \item{object}{a \code{riley} object.}
  \item{\dots}{arguments to be passed on to other functions}
}

\details{The variance-covariance matrix is obtained from the inverse Hessian as provided by \code{optim}.}

\value{A matrix of the estimated covariances between the parameter estimates in the Riley model: logit of sensitivity (mu1), logit of false positive rate (\code{mu2}), additional variation of \code{mu1} beyond sampling error (\code{psi1}), additional variation of \code{mu2} beyond sampling error (\code{psi2}) and a transformation of the correlation between \code{psi1} and \code{psi2} (\code{rhoT}). The original correlation is given as \code{inv.logit(rhoT)*2-1}.}

\note{A warning message is casted when the Hessian matrix contains negative eigenvalues. This implies that the identified minimum for the (restricted) negative log-likelihood is a saddle point, and that the solution is therefore not optimal.}

\references{
Riley, RD., Thompson, JR., & Abrams, KR. (2008). \dQuote{An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown.} \emph{Biostatistics}, \bold{9}, 172--186.
}

\author{Thomas Debray <thomas.debray@gmail.com>}



\seealso{
\code{\link{riley}}
}

\keyword{riley}
\keyword{bivariate}
