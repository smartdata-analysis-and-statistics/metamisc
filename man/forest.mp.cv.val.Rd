% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metapred_measures.R
\name{forest.mp.cv.val}
\alias{forest.mp.cv.val}
\alias{forest.perf}
\title{Forest plot of a validation object.}
\usage{
\method{forest}{mp.cv.val}(x, perfFUN = 1, method = "REML", xlab = NULL, ...)
}
\arguments{
\item{x}{An \code{mp.cv.val} or \code{perf} object.}

\item{perfFUN}{Numeric or character. Which performance statistic should be plotted? Defaults to the first.}

\item{method}{character string specifying whether a fixed- or a random-effects model should be used to summarize the
prediction model performance. A fixed-effects model is fitted when using method="FE". Random-effects models are fitted 
by setting method equal to one of the following: "DL", "HE", "SJ", "ML", "REML", "EB", "HS", or "GENQ". Default is "REML".}

\item{xlab}{Label on x-axis. Defaults to the name of the performance function.}

\item{...}{Other arguments passed to plotting internals. E.g. \code{title}. See \link{forest.default} for details.}
}
\description{
Draw a forest plot of the performance of an internally-externally cross-validated model.
}
\author{
Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>
}
