context("valmeta 3. bayesian meta-analysis of c-statistic")
skip_on_cran()


library(metafor)
library(dplyr)
data(EuroSCORE)

test_that("Bayesian random effect meta-analysis of concordance statistic", {
  fit <- valmeta(cstat = c.index, cstat.se = se.c.index, cstat.cilb = c.index.95CIl,
                 cstat.ciub = c.index.95CIu, cstat.cilv = 0.95, N = n, O = n.events, 
                 data = EuroSCORE, method = "BAYES", slab = Study)
  expect_equal(as.character(class(fit$fit)), "runjags")
  
})


test_that("Bayesian random effect meta-analysis of concordance statistic using normality/identity model", {
  
  pars <- list(model.cstat = "normal/identity")
  
  fit <- valmeta(cstat = c.index, cstat.se = se.c.index, cstat.cilb = c.index.95CIl,
                 cstat.ciub = c.index.95CIu, cstat.cilv = 0.95, N = n, O = n.events, 
                 data = EuroSCORE, method = "BAYES", slab = Study, pars = pars)
  
  expect_equal(as.character(class(fit$fit)), "runjags")
  
})