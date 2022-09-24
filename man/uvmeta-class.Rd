\name{uvmeta-class}
\alias{uvmeta-clas}

%\alias{slice,mle2-method}
%\alias{summary,mle2-method}
%\alias{update,mle2-method}
%\alias{vcov,mle2-method}
%\alias{deviance,mle2-method}
%\alias{coerce,mle,mle2-method}
%\alias{formula,mle2-method}
%\alias{stdEr}
%\alias{stdEr,mle2-method}
\title{Class "uvmeta". Result of a univariate meta-analysis.}
\description{This class encapsulates results of a univariate meta-analysis.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{\link{uvmeta}}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{(language) The call to \code{\link{uvmeta}}.}
    \item{\code{data}:}{(data frame) The data used for the meta-analysis. }
    \item{\code{results}:}{(data frame) Contains the pooled effect size (\code{mu}), the between-study variability (\code{tausq}), Cochran's Q statistic (\code{Q}) and Higgins' and Thompson's I square statistic (\code{Isq}). For each estimate, error variances are provided with predefined confidence (\code{method="MOM"}) or credibility (\code{method="bayes"}) intervals. }
    \item{\code{model}:}{(character) The meta-analysis model used.}
    \item{\code{method}:}{(character) The estimator used.}
    \item{\code{na.action}:}{(character) Information from the action which was applied to object if NAs were handled specially, or NULL.}
    \item{\code{df}:}{(numeric) Degrees of freedom.}
    \item{\code{numstudies}:}{(numeric) The amount of studies used in the meta-analysis.}
    \item{\code{pred.int}:}{(data frame) A prediction interval, predicting in what range future effect sizes will fall given what has already been observed (based on a Student's t-distribution, cfr. Riley 2011)}
    \item{\code{formula}:}{(character) If a formula was specified, a
      character vector giving the formula and parameter specifications.}
}
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(object = "uvmeta")}: Print object summary.}
    \item{forest}{\code{signature(object = "uvmeta")}: Plot a forest plot with the summary estimate.}
    \item{summary}{\code{signature(object = "uvmeta")}: Generate object summary.}

  }
}

\examples{
data(Collins)

#Extract effect size and error variance
r <- Collins$logOR
vars <- Collins$SE**2

#Frequentist random-effects meta-analysis
fit1 <- uvmeta(r,vars)

#Extract results
fit1$results
}
\keyword{classes}