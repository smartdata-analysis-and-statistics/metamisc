% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metapred.R
\name{metapred}
\alias{metapred}
\title{Generalized Stepwise Regression for Prediction Models in Clustered Data}
\usage{
metapred(
  data,
  strata,
  formula,
  estFUN = "glm",
  scope = NULL,
  retest = FALSE,
  max.steps = 1000,
  center = FALSE,
  recal.int = FALSE,
  cvFUN = NULL,
  cv.k = NULL,
  metaFUN = NULL,
  meta.method = NULL,
  predFUN = NULL,
  perfFUN = NULL,
  genFUN = NULL,
  selFUN = "which.min",
  gen.of.perf = "first",
  ...
)
}
\arguments{
\item{data}{data.frame containing the data. Note that \code{metapred} removes observations with missing data \emph{listwise}
for all variables in \code{formula} and \code{scope}, to ensure that the same data is used in each model in each step. The
outcome variable should be \code{numeric} or coercible to numeric by as.numeric().}

\item{strata}{Character to specify the name of the strata (e.g. studies or clusters) variable}

\item{formula}{\code{formula} of the first model to be evaluated. \code{metapred} will start at \code{formula} and update it
using terms of \code{scope}. Defaults to full main effects model, where the first column in \code{data} is assumed to be
the outcome and all remaining columns (except \code{strata}) predictors. See \link[stats]{formula} for formulas in general.}

\item{estFUN}{Function for estimating the model in the first stage. Currently "lm", "glm" and "logistfirth" are supported.}

\item{scope}{\code{formula}. The difference between \code{formula} and \code{scope} defines the range of models examined in the 
stepwise search. Defaults to NULL, which leads to the intercept-only model. If \code{scope} is not nested in \code{formula}, 
this implies backwards selection will be applied (default). If \code{scope} is nested in \code{formula}, this implies forward 
selection will be applied. If equal, no stepwise selection is applied.}

\item{retest}{Logical. Should added (removed) terms be retested for removal (addition)? \code{TRUE} implies bi-directional 
stepwise search.}

\item{max.steps}{Integer. Maximum number of steps (additions or removals of terms) to take. Defaults to 1000, which is
essentially as many as it takes. 0 implies no stepwise selection.}

\item{center}{logical. Should numeric predictors be centered around the cluster mean?}

\item{recal.int}{Logical. Should the intercept be recalibrated in each validation?}

\item{cvFUN}{Cross-validation method, on the study (i.e. cluster or stratum) level. "l1o" for leave-one-out cross-validation 
(default). "bootstrap" for bootstrap. Or "fixed", for one or more data sets which are only used for validation. A user written 
function may be supplied as well.}

\item{cv.k}{Parameter for cvFUN. For \code{cvFUN="bootstrap"}, this is the number of bootstraps. For \code{cvFUN="fixed"}, 
this is a vector of the indices of the (sorted) data sets. Not used for \code{cvFUN="l1o"}.}

\item{metaFUN}{Function for computing the meta-analytic coefficient estimates in two-stage MA. 
By default, \link[metafor]{rma.uni}, from the metafor package is used. Default settings are univariate random effects, 
estimated with "REML". Method can be passed trough the \code{meta.method} argument.}

\item{meta.method}{Name of method for meta-analysis. Default is "REML". For more options see \link[metafor]{rma.uni}.}

\item{predFUN}{Function for predicting new values. Defaults to the predicted probability of the outcome, using the link 
function of \code{glm()} or \code{lm()}.}

\item{perfFUN}{Function for computing the performance of the prediction models. Default: mean squared error 
(\code{perfFUN="mse"}, aka Brier score for binomial outcomes).Other options are \code{"var.e"} (variance of prediction error), 
\code{"auc"} (area under the curve),
\code{"cal_int"} (calibration intercept), and \code{"cal_slope"} (multiplicative calibration slope) and \code{"cal_add_slope"}
(additive calibration slope), or a \code{list} of these, where only the first is used for model selection.}

\item{genFUN}{Function or \code{list} of named functions for computing generalizability of the performance. 
Default: \code{rema}, summary statistic of a random effects meta-analysis. Choose \code{"rema_tau"} for heterogeneity 
estimate of a random effects meta-analysis,  \code{genFUN="abs_mean"} for (absolute) mean,  
\code{coefficient_of_variation} for the coefficient of variation. If a \code{list} containing these, only the first is used 
for model selection.}

\item{selFUN}{Function for selecting the best method. Default: lowest value for \code{genFUN}. Should be set to
"which.max" if high values for \code{genFUN} indicate a good model.}

\item{gen.of.perf}{For which performance measures should generalizability measures be computed? \code{"first"} (default) for 
only the first. \code{"respective"} for matching the generalizability measure to the performance measure on the same location
in the list. \code{"factorial"} for applying all generalizability measures to all performance estimates.}

\item{...}{To pass arguments to estFUN (e.g. family = "binomial"), or to other FUNctions.}
}
\value{
A list of class \code{metapred}, containing the final model in \code{global.model}, and the stepwise
tree of estimates of the coefficients, performance measures, generalizability measures in \code{stepwise}.
}
\description{
Generalized stepwise regression for obtaining a prediction model that is validated 
with (stepwise) internal-external cross-validation, in or to obtain adequate 
performance across data sets. Requires data from individuals in multiple studies.
}
\details{
Use \link{subset.metapred} to obtain an individual prediction model from a \code{metapred} object.

 Note that \code{formula.changes} is currently unordered; it does not represent the order of changes in the stepwise 
 procedure.
 
 \code{metapred} is still under development, use with care.
}
\examples{
data(DVTipd)

\dontrun{
# Explore heterogeneity in intercept and assocation of 'ddimdich'
glmer(dvt ~ 0 + cluster + (ddimdich|study), family = binomial(), data = DVTipd)
}

# Scope
f <- dvt ~ histdvt + ddimdich + sex + notraum

# Internal-external cross-validation of a pre-specified model 'f'
fit <- metapred(DVTipd, strata = "study", formula = f, scope = f, family = binomial)
fit

# Let's try to simplify model 'f' in order to improve its external validity
metapred(DVTipd, strata = "study", formula = f, family = binomial)

# We can also try to build a generalizable model from scratch

\dontrun{
# Some additional examples:
metapred(DVTipd, strata = "study", formula = dvt ~ 1, scope = f, family = binomial) # Forwards
metapred(DVTipd, strata = "study", formula = f, scope = f, family = binomial) # no selection
metapred(DVTipd, strata = "study", formula = f, max.steps = 0, family = binomial) # no selection
metapred(DVTipd, strata = "study", formula = f, recal.int = TRUE, family = binomial)
metapred(DVTipd, strata = "study", formula = f, meta.method = "REML", family = binomial)
}
# By default, metapred assumes the first column is the outcome.
newdat <- data.frame(dvt=0, histdvt=0, ddimdich=0, sex=1, notraum=0)
fitted <- predict(fit, newdata = newdat)

}
\references{
Debray TPA, Moons KGM, Ahmed I, Koffijberg H, Riley RD. A framework for developing, implementing, 
and evaluating clinical prediction models in an individual participant data meta-analysis. 
\emph{Stat Med}. 2013;32(18):3158-80.

de Jong VMT, Moons KGM, Eijkemans MJC, Riley RD, Debray TPA. Developing more generalizable prediction models from pooled 
studies and large clustered data sets. \emph{Stat Med}. 2021;40(15):3533--59.

Riley RD, Tierney JF, Stewart LA. Individual participant data meta-analysis: a handbook for healthcare research. 
Hoboken, NJ: Wiley; 2021. ISBN: 978-1-119-33372-2.

Schmid CH, Stijnen T, White IR. Handbook of meta-analysis. First edition. Boca Raton: Taylor and Francis; 2020. ISBN: 978-1-315-11940-3.
}
\seealso{
\code{\link{forest.metapred}}  for generating a forest plot of prediction model performance
}
\author{
Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>
}
