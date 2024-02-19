\name{Tzoulaki}
\alias{Tzoulaki}
\docType{data}
\title{
The incremental value of cardiovascular risk factors
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Tzoulaki et al. (2009) reviewed studies that evaluated various candidate prognostic factors in their ability to improve prediction of coronary heart disease (CHD) or other outcomes beyond what the Framingham risk score (FRS) can achieve.

}
\usage{data("Tzoulaki")}
\format{
  A data frame containing data from 27 studies on the following 2 variables.
  \describe{
    \item{\code{PubmedID}}{a character vector with the Pubmed ID of the study}
    \item{\code{N}}{a numeric vector describing the study size}
    \item{\code{N.events}}{a numeric vector describing the observed number of events}
    \item{\code{FRS.orig.refitted}}{a boolean vector describing whether the coefficients of the original Framingham Risk Score (FRS) were re-estimated}
    \item{\code{FRS.modif.refitted}}{a boolean vector describing whether the coefficients of the modified Framingham Risk Score were re-estimated}
    \item{\code{predictors}}{a character vector indicating which new risk factor(s) were included in the modified FRS }
    \item{\code{outcome}}{a character vector indicating the primary outcome being predicted}
    \item{\code{AUC.orig}}{a numeric vector describing the Area under the ROC curve (AUC) of the original FRS model}
    \item{\code{AUC.orig.CIl}}{a numeric vector describing the lower boundary of the 95\% confidence interval of the AUC of the original FRS model}
    \item{\code{AUC.orig.CIu}}{a numeric vector describing the upper boundary of the 95\% confidence interval of the AUC of the original FRS model}
    \item{\code{AUC.modif}}{a numeric vector describing the Area under the ROC curve (AUC) of the modified FRS model that includes one or more new risk factors}
    \item{\code{AUC.modif.CIl}}{a numeric vector describing the lower boundary of the 95\% confidence interval of the AUC of the modified FRS model}
    \item{\code{AUC.modif.CIu}}{a numeric vector describing the upper boundary of the 95\% confidence interval of the AUC of the modified FRS model}
    \item{\code{pval.AUCdiff}}{a numeric vector with the p-value of the difference between \code{AUC.orig} and \code{AUC.modif}}
    \item{\code{sign.AUCdiff}}{a boolean vector indicating whether the difference between \code{AUC.orig} and \code{AUC.modif} is below 0.05}
  }
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Tzoulaki I, Liberopoulos G, Ioannidis JPA. Assessment of claims of improved prediction beyond the Framingham risk score.  \emph{JAMA}. 2009 Dec 2;302(21):2345–52. 

}
\examples{
data(Tzoulaki)
## maybe str(Tzoulaki) ; plot(Tzoulaki) ...
}
\keyword{datasets}
