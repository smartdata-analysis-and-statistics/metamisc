\name{Zhang}
\alias{Zhang}
\docType{data}
\title{
Meta-analysis of the prognostic role of hormone receptors in endometrial cancer
}
\description{
This dataset comprises the results from 16 studies assessing the prognostic role of human epidermal growth factor receptor 2 (HER2) in endometrial cancer. These studies were previously identified in a systematic review by Zhang et al. to evaluate the overall risk of several hormone receptors for endometrial cancer survival. 
}
\usage{data("Zhang")}
\format{
  A data frame with 20 observations on the following 10 variables.
  \describe{
    \item{\code{Study}}{a factor with 16 levels to indicate the study}
    \item{\code{PrimaryAuthor}}{a factor indicating the first author's last name}
    \item{\code{year}}{a numeric vector indicating the publication year}
    \item{\code{Country}}{a factor indicating the source country of the study data}
    \item{\code{Disease}}{a factor indicating the studied disease. Possible levels are \code{EC} (endometrial cancer), \code{EEC} (endometrioid endometrial cancer) and \code{UPSC} (uterine papillary serous carcinoma)}
    \item{\code{N}}{a numeric vector describing the total sample size of each study}
    \item{\code{HR}}{a numeric vector describing the estimated hazard ratio of each study}
    \item{\code{HR.025}}{a numeric vector describing the lower boundary of the 95\% confidence interval of \code{HR}}
    \item{\code{HR.975}}{a numeric vector describing the upper boundary of the 95\% confidence interval of \code{HR}}
    \item{\code{outcome}}{a factor indicating the studied outcome. Possible levels are \code{OS} (overall survival) and \code{PFS} (progression-free survival)}
  }
}
\details{
Eligible studies were identified by searching the PubMed and EMBASE databases for publications from 1979 to May 2014. Data were collected from studies comparing overall survival or progression-free survival in patients with elevated levels of human epidermal growth factor receptor 2 with those in patients with lower levels. 
}
\source{
Zhang Y, Zhao D, Gong C, Zhang F, He J, Zhang W, et al. Prognostic role of hormone receptors in endometrial cancer: a systematic review and meta-analysis. World J Surg Oncol. 2015 Jun 25;13:208. 
}
\references{
Riley RD, Jackson D, Salanti G, Burke DL, Price M, Kirkham J, et al. Multivariate and network meta-analysis of multiple outcomes and multiple treatments: rationale, concepts, and examples. BMJ. 2017 13;358:j3932. 
}
\examples{
data(Zhang)

# Display the hazard ratios for overall survival in a forest plot
ds <- subset(Zhang, outcome=="OS")
with(ds, forest(theta = HR, theta.ci.lb = HR.025, theta.ci.ub = HR.975, 
                theta.slab = Study, xlab = "Hazard ratio of HER2 versus OS", refline = 1))
}
\keyword{datasets}
