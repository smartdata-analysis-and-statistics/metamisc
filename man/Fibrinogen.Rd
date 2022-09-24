\name{Fibrinogen}
\alias{Fibrinogen}
\docType{data}
\title{
Meta-analysis of the association between plasma fibrinogen concentration and the risk of coronary heath disease
}
\description{
The \code{Fibrinogen} data set is a meta-analysis of 31 studies in which the association between plasma fibrinogen concentration and the risk of coronary heath disease (CHD) was estimated.
}
\usage{data("Fibrinogen")}
\format{
  A data frame with 5 variables:
  \describe{
    \item{\code{N.total}}{a numeric vector describing the total number of patients for each study}
    \item{\code{N.events}}{a numeric vector describing the number of observed events within each study}
    \item{\code{HR}}{a numeric vector describing the estimated hazard ratio of each study}
    \item{\code{HR.025}}{a numeric vector describing the lower boundary of the 95\% confidence interval of \code{HR}}
    \item{\code{HR.975}}{a numeric vector describing the upper boundary of the 95\% confidence interval of \code{HR}}
  }
}

\source{
Fibrinogen Studies Collaboration. Collaborative meta-analysis of prospective studies of plasma fibrinogen and cardiovascular disease. \emph{Eur J Cardiovasc Prev Rehabil}. 2004 Feb;11(1):9-17. 

Thompson S, Kaptoge S, White I, Wood A, Perry P, Danesh J, et al. Statistical methods for the time-to-event analysis of individual participant data from multiple epidemiological studies. \emph{Int J Epidemiol}. 2010 Oct;39(5):1345-59. 
}

\examples{
data(Fibrinogen)
## maybe str(Fibrinogen) ; plot(Fibrinogen) ...
}
\keyword{datasets}
