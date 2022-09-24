context("valmeta 4. frequentist meta-analysis of total O:E ratio")
skip_on_cran()

library(lme4)
library(dplyr)
data(EuroSCORE)

test_that("Test whether Using the data argument yields the same results", {
  
  # valmeta using O and E
  fit1 <- valmeta(measure="OE", O=n.events, E=e.events, data=EuroSCORE)     
  fit2 <- with(EuroSCORE, valmeta(measure="OE", O=n.events, E=e.events))
  expect_equal(fit1$numstudies, fit2$numstudies)
  expect_equal(fit1$est, fit2$est)
  expect_equal(fit1$ci.lb, fit2$ci.lb)
  expect_equal(fit1$ci.ub, fit2$ci.ub)
  expect_equal(fit1$pi.lb, fit2$pi.lb)
  expect_equal(fit1$pi.ub, fit2$pi.ub)
  
  # valmeta using O, E and N
  fit3 <- valmeta(measure="OE", O=n.events, E=e.events, N=n,data=EuroSCORE)
  fit4 <- with(EuroSCORE, valmeta(measure="OE", O=n.events, E=e.events, N=n))
  expect_equal(fit3$numstudies, fit4$numstudies)
  expect_equal(fit3$est, fit4$est)
  expect_equal(fit3$ci.lb, fit4$ci.lb)
  expect_equal(fit3$ci.ub, fit4$ci.ub)
  expect_equal(fit3$pi.lb, fit4$pi.lb)
  expect_equal(fit3$pi.ub, fit4$pi.ub)
  
  # valmeta using Po, Pe and N
  fit5 <- valmeta(measure="OE", Po=Po, Pe=Pe, N=n, data=EuroSCORE)
  fit6 <- with(EuroSCORE, valmeta(measure="OE", Po=Po, Pe=Pe, N=n))
  expect_equal(fit5$numstudies, fit6$numstudies)
  expect_equal(fit5$est, fit6$est)
  expect_equal(fit5$ci.lb, fit6$ci.lb)
  expect_equal(fit5$ci.ub, fit6$ci.ub)
  expect_equal(fit5$pi.lb, fit6$pi.lb)
  expect_equal(fit5$pi.ub, fit6$pi.ub)
  
})

test_that("Test whether changing the parameters works", {
  
  fit.model <- "poisson/log"
  fit <- with(EuroSCORE, valmeta(measure="OE", O=n.events, E=e.events, method="ML", test="z",
                                 pars=list(model.oe=fit.model)))
  
  expect_equal(fit$model, fit.model)
  
  fit.model <- "normal/identity"
  fit <- with(EuroSCORE, valmeta(measure="OE", O=n.events, E=e.events,
                                 pars=list(model.oe=fit.model)))
  
  expect_equal(fit$model, fit.model)
})

test_that("Fixed effect meta-analysis of total O:E ratio works (Poisson distribution)", {

  # Let's ignore clustering of studies in this test

  fit1.glm <- glm (n.events ~ 1, offset = log(e.events),
                   family = poisson(link = "log"), data = EuroSCORE)

  fit1.valmeta <- with(EuroSCORE, valmeta(measure="OE",
                                          O=n.events,
                                          E=e.events,
                                          test = "z",
                                          method="FE", slab = Study,
                                          pars=list(model.oe="poisson/log")))

  # Now try the same but omit study labels
  fit2.valmeta <- with(EuroSCORE, valmeta(measure="OE",
                                          O=n.events,
                                          E=e.events,
                                          test = "z",
                                          method="FE",
                                          pars=list(model.oe="poisson/log")))

  expect_equal(fit1.valmeta$est, fit2.valmeta$est, exp(as.numeric(coefficients(fit1.glm))))
})

test_that("Random effects meta-analysis of total O:E ratio using Poisson distribution", {

  # Let's ignore clustering of studies in this test
  ds <- EuroSCORE
  ds$Study <- c(1:dim(ds)[1]) #Re-assign study labels

  fit1.lme4 <- glmer(n.events ~ 1 | Study, offset = log(e.events),
                     family = poisson(link = "log"), data = ds)

  fit1.valmeta <- with(EuroSCORE, valmeta(measure = "OE",
                                          O = n.events,
                                          E = e.events,
                                          test = "z",
                                          method = "ML", 
                                          slab = Study,
                                          pars = list(model.oe="poisson/log")))

  # Now try the same but omit study labels
  fit2.valmeta <- with(EuroSCORE, valmeta(measure="OE",
                                          O=n.events,
                                          E=e.events,
                                          test = "z",
                                          method="ML",
                                          pars=list(model.oe="poisson/log")))

  expect_equal(fit1.valmeta$est, fit2.valmeta$est, exp(as.numeric(fixef(fit1.lme4))))

  # valmeta should still work if some studies cannot be used
  O <- EuroSCORE$n.events
  O[1:10] <- NA
  fit3.valmeta <- valmeta(measure="OE",
                          O=O,
                          E=EuroSCORE$e.events,
                          slab = EuroSCORE$Study,
                          test = "z", method = "ML",
                          pars=list(model.oe="poisson/log"))
  expect_equal(fit3.valmeta$numstudies, sum(!is.na(O)))

  # valmeta should still work if observed event counts are not properly rounded
  O <- EuroSCORE$n.events
  O[1:10] <- O[1:10] + 0.2
  fit4.valmeta <- valmeta(measure="OE",
                          O=O,
                          E=EuroSCORE$e.events,
                          slab = EuroSCORE$Study,
                          test = "z", method = "ML",
                          pars=list(model.oe="poisson/log"))


})

test_that("Generating a plot", {

  expect_is(fit1.valmeta <- valmeta(measure = "OE",
                                          O = n.events,
                                          E = e.events,
                                          test = "z",
                                          method="FE", 
                                    slab = Study, 
                                    data = EuroSCORE)
            , "valmeta")

  plot(fit1.valmeta)


})

test_that("Test storage of results", {
  
  # Valmeta using cstat, se.cstat, N, and O
  fit1 <- valmeta(measure = "OE", N=n, O=n.events, E =  e.events,
                 slab=Study, data=EuroSCORE,
                 method     = "ML",
                 pars       = list(model.oe = "poisson/log"))
  
  expect_true("data.frame" %in% class(fit1$data))
  expect_true(all(c("theta", "theta.se", "theta.cilb", "theta.ciub") %in% colnames(fit1$data)))
  
})




