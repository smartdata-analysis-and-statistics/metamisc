\name{Daniels}
\alias{Daniels}
\docType{data}

\title{Daniels and Hughes data}

\description{Data frame with treatment differences in CD4 cell count.}

\usage{data("Daniels")}


\format{A data frame with 15 observations on the following 2 variables.
  \describe{
  \item{\code{Y1}}{Treatment differences for the log hazard ratio for the development of AIDS or death over 2 years.}  
  \item{\code{vars1}}{Error variances of \code{Y1}.}  
  \item{\code{Y2}}{Difference in mean change in CD4 cell count between baseline and 6 month for studies of the AIDS Clinical Trial Group}  
  \item{\code{vars2}}{Error variances of \code{Y2}.}  
  }
}

\details{The \code{Daniels} data comprises 15 phase II/III randomized clinical trials of the HIV Disease Section of the Adult AIDS Clinical Trials Group of the National Institutes of Health, which had data available as of May 1996, which had at least six months of follow-up on some patients and in which at least one patient developed AIDS or died. The data were previously used
by Daniels and Hughes (1997) to assess whether the change in CD4 cell count is a surrogate for time to either development of AIDS or death in drug trials of patients with HIV. 
}

\source{
Daniels MJ, Hughes MD. Meta-analysis for the evaluation of potential surrogate markers. \emph{Statistics in Medicine} 1997; \bold{16}: 1965--1982.
}

\keyword{datasets}
