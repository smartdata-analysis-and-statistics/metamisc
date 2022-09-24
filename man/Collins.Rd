\name{Collins}
\alias{Collins}
\docType{data}

\title{Collins data}

\description{A meta-analysis of nine clinical trials investigating the effect of taking diuretics during pregnancy on the risk of pre-eclampsia. 
}
\usage{data(Collins)}
\format{
  A data frame with 9 observations on the following 2 variables.
  \describe{
    \item{\code{logOR}}{a numeric vector with treatment effect sizes (log odds ratio)}
    \item{\code{SE}}{a numeric vector with the standard error of the treatment effect sizes}
  }
}

\source{
Collins, R., Yusuf, S., Peto, R. Overview of randomised trials of diuretics in pregnancy. \emph{British Medical Journal} 1985, \bold{290}, 17--23.

Hardy, R.J. Thompson, S.G. A likelihood approach to meta-analysis with random effects. \emph{Statistics in Medicine} 1996; \bold{15}:619--629.
}

\examples{
data(Collins)
}
\keyword{datasets}
