\name{impact}
\alias{impact}
\docType{data}
\title{IMPACT data}
\description{
The IMPACT dataset comprises 15 studies of patients suffering from traumatic brain injury, including individual patient data from 11 randomized controlled trials and four observational studies. 
}
\usage{data("impact")}
\format{
  A data frame with 11022 observations on the following 11 variables.
  \describe{
    \item{\code{name}}{Name of the study}
    \item{\code{type}}{Type of study, RCT: randomized controlled trial,OBS: observational cohort}
    \item{\code{age}}{Age of the patient}
    \item{\code{motor_score}}{Glasgow Coma Scale motor score}
    \item{\code{pupil}}{Pupillary reactivity}
    \item{\code{ct}}{Marshall Computerized Tomography classification}
    \item{\code{hypox}}{Hypoxia (0=no, 1=yes)}
    \item{\code{hypots}}{Hypotension (0=no, 1=yes)}
    \item{\code{tsah}}{Traumatic subarachnoid hemorrhage (0=no, 1=yes)}
    \item{\code{edh}}{Epidural hematoma (0=no, 1=yes)}
    \item{\code{mort}}{6-month mortality (0=alive, 1=dead)}
  }
}
\details{
The included studies were part of the IMPACT project, where a total of 25 prognostic factors were considered for prediction of 6-month mortality. Missing values were imputed using the study as a fixed effect in the imputation model (Steyerberg et al, 2008).
}
\source{
Steyerberg EW, Nieboer D, Debray TPA, Van Houwelingen JC. Assessment of heterogeneity in an individual participant data meta-analysis of prediction models: An overview and illustration. \emph{Stat Med}. 2019;38(22):4290--309. 
}
\references{
Murray GD, Butcher I, McHugh GS, et al. Multivariable prognostic analysis in traumatic brain injury: results from the IMPACT study. \emph{J Neurotrauma}. 2007;24(2):329--337.

Steyerberg EW, Mushkudiani N, Perel P, et al. Predicting outcome after traumatic brain injury: development and international validation of prognostic scores based on admission characteristics. \emph{PLOS Med}. 2008;5(8):e165.
}
\examples{
data(impact)

by(impact, impact$name, summary)

# Plot the distribution of age by study
library(ggplot2)
e <- ggplot(impact, aes(x = name, y = age))
e + geom_violin(aes(fill = type), trim = FALSE)  + 
    theme(axis.text.x = element_text(angle = 45)) + 
    xlab("Study")
}
\keyword{datasets}
