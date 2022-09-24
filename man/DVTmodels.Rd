\name{DVTmodels}
\alias{DVTmodels}
\docType{data}
\title{Risk prediction models for diagnosing Deep Venous Thrombosis (DVT)
}
\description{Previously published prediction models for predicting the presence of DVT. }

\usage{data(DVTmodels)}

\format{
An object of the class \code{litmodels} with the following information for each literature model: the study-level descriptives ("\code{descriptives}"), the regression coefficient or weight for each predictor ("\code{weights}") and the error variance for each regression coefficient or weight ("\code{weights.var}"). 
}

\details{
Previously, several models (Gagne, Oudega) and score charts (Wells, modified Wells, and Hamilton) have been published for evaluating the presence of DVT in suspected patients. These models combine information on mulitple predictors into a weighted sum, that can subsequently be used to obtain estimates of absolute risk. See \code{DVTipd} for more information on the predictors.
}
\source{
Wells PS, Anderson DR, Bormanis J, Guy F, Mitchell M, Gray L, Clement C, Robinson KS, Lewandowski B. Value of assessment of pretest probability of deep-vein thrombosis in clinical management. \emph{Lancet} 1997; \bold{350}(9094):1795--1798. DOI: 10.1016/S0140-6736(97)08140-3.

Wells PS, Anderson DR, Rodger M, Forgie M, Kearon C, Dreyer J, Kovacs G, Mitchell M, Lewandowski B, Kovacs MJ. Evaluation of D-dimer in the diagnosis of suspected deep-vein thrombosis. \emph{The New England Journal of Medicine} 2003; \bold{349}(13):1227--1235. DOI: 10.1056/NEJMoa023153.

Gagne P, Simon L, Le Pape F, Bressollette L, Mottier D, Le Gal G. Clinical prediction rule for diagnosing deep vein thrombosis in primary care. \emph{La Presse Medicale} 2009; \bold{38}(4):525--533. DOI: 10.1016/j.lpm.2008.09.022.

Subramaniam RM, Snyder B, Heath R, Tawse F, Sleigh J. Diagnosis of lower limb deep venous thrombosis in emergency department patients: performance of Hamilton and modified Wells scores. \emph{Annals of Emergency Medicine} 2006; \bold{48}(6):678--685. DOI: 10.1016/j.annemergmed.2006.04.010.

Oudega R, Moons KGM, Hoes AW. Ruling out deep venous thrombosis in primary care. a simple diagnostic algorithm including D-dimer testing. \emph{Thrombosis and Haemostasis} 2005; \bold{94}(1):200--205. DOI: 10.1160/TH04-12-0829.
}

\references{
Debray TPA, Koffijberg H, Nieboer D, Vergouwe Y, Steyerberg EW, Moons KGM. Meta-analysis and aggregation of multiple published prediction models. \emph{Stat Med}. 2014 Jun 30;33(14):2341--62. 

Debray TPA, Koffijberg H, Vergouwe Y, Moons KGM, Steyerberg EW. Aggregating published prediction models with individual participant data: a comparison of different approaches. \emph{Stat Med}. 2012 Oct 15;31(23):2697--712. 

}

\seealso{
\code{\link{DVTipd}}
}
\examples{
data(DVTmodels)
}
\keyword{datasets}
