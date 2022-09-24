\name{inv.logit}
\alias{inv.logit}

\title{
  Apply the inverse logit tranformation
}
\description{
  Transforms a linear predictor into a probability.
}

\usage{
  inv.logit(x)
}

\arguments{
  \item{x}{A vector of numerics (between -Inf and Inf)}
}

\value{
  A vector of numerics between 0 and 1.
}
\author{
  Thomas Debray <thomas.debray@gmail.com>
}

\seealso{
  \code{\link{logit}}
}

\keyword{logit}
\keyword{inverse}
\keyword{transformation}
