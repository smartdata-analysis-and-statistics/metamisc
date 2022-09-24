\name{Framingham}
\alias{Framingham}
\docType{data}
\title{
Predictive performance of the Framingham Risk Score in male populations
}
\description{
This data set contains estimates on the performance of the Framingham model for predicting coronary heart disease in male populations (Wilson 1998). Results are based on the original development study and 20 validations identified by Damen \emph{et al} (BMC Med, 2017).
}
\usage{data("Framingham")}
\format{
  A data frame with 24 observations on the following 19 variables.
  \describe{
    \item{\code{AuthorYear}}{a vector describing the study authors}
    \item{\code{n}}{a numeric vector with the total number of patients on which performance estimates are based}
    \item{\code{n.events}}{a numeric vector with the total number of observed events}
    \item{\code{c.index}}{a numeric vector with the estimated concordance statistic of each validation}
    \item{\code{se.c.index}}{a numeric vector with the standard error of the concordance statistics}
    \item{\code{c.index.95CIl}}{a numeric vector with the lower bound of the 95\% confidence interval of the estimated concordance statistics}
    \item{\code{c.index.95CIu}}{a numeric vector with the upper bound of the 95\% confidence interval of the estimated concordance statistics}
    \item{\code{Po}}{a numeric vector with the overall observed event probability of each validation}
    \item{\code{Pe}}{a numeric vector with the overall expected event probability of each validation}
    \item{\code{t.val}}{a numeric vector describing the time period in which predictive performance was assessed for each validation}
    \item{\code{mean_age}}{a numeric vector describing the mean age of the patients}
    \item{\code{sd_age}}{a numeric vector with the spread of the age of the patients}
    \item{\code{mean_SBP}}{a numeric vector with the mean systolic blood pressure in the validation studies (mm Hg)}
    \item{\code{sd_SBP}}{a numeric vector with the spread of systolic blood pressure in the validation studies}
    \item{\code{mean_total_cholesterol}}{a numeric vector with the mean total cholesterol in the validation studies (mg/dL)}
    \item{\code{sd_total_cholesterol}}{a numeric vector with the spread of total cholesterol in the validation studies}
    \item{\code{mean_hdl_cholesterol}}{a numeric vector with the mean high-density lipoprotein cholesterol in the validation studies (mg/dL)}
    \item{\code{sd_hdl_cholesterol}}{a numeric vector with the spread of high-density lipoprotein cholesterol in the validation studies}
    \item{\code{pct_smoker}}{a numeric vector with the percentage smokers in the validation studies}
  }
}
\details{
The Framingham Risk Score allows physicians to predict 10-year coronary heart disease (CHD) risk in patients without overt CHD. It was developed in 1998 from a middle-aged white population sample, and has subsequently been validated across different populations. The current dataset contains the original (internal validation) results, as well as 23 external validations which were identified through a systematic review. In this review, studies were eligible for inclusion if they described the validation of the original Framingham model and assessed its performance for fatal or nonfatal CHD in males from a general population setting.
}
\source{
Damen JAAG, Hooft L, Schuit E, Debray TPA, Collins GS, Tzoulaki I, et al. Prediction models for cardiovascular disease risk in the general population: systematic review. \emph{BMJ}. 2016;i2416. 

Damen JAAG, Pajouheshnia R, Heus P, Moons KGM, Reitsma JB, Scholten RJPM, et al. Performance of the Framingham risk models and Pooled Cohort Equations: a systematic review and meta-analysis. \emph{BMC Med}. 2017;17(1):109. 

Wilson PW, D'Agostino RB, Levy D, Belanger AM, Silbershatz H, Kannel WB. Prediction of coronary heart disease using risk factor categories. \emph{Circulation}. 1998; \bold{97}(18):1837--47. 

}
\examples{
data(Framingham)
}
\keyword{datasets}
