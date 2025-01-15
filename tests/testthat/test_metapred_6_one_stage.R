context("metapred 6. One-stage.")

set.seed(7318)
n <- 100
y <- c(rep(0, 2 * n), rep(1, 2 * n))
x <- rnorm(length(y), y, 1)
y <- factor(y)
k <- rep(1:8, 5* n)
d <- data.frame(y, x, k)

test_that("metapred can estimate one-stage fixed effect models", { 
  skip_on_cran()
  f <- y ~ x
  mp_fe <- metapred(d, "k", formula = f, scope = f, family = binomial, estFUN = glm, two.stage = F)
  expect_is(mp_fe, "metapred")
})


test_that("metapred one-stage can extract fitted values.", {
  skip_on_cran()
  
  f <- y ~ x
  mp_fe <- metapred(d, "k", formula = f, scope = f, family = binomial, estFUN = glm, two.stage = F)
  
  fitted_values_stratified <- fitted(mp_fe)
  expect_vector(fitted_values_stratified)
  expect_length(fitted_values_stratified, length(unique(k)))
  
  fitted_values_unlisted <- fitted(mp_fe, as.stratified = FALSE)
  expect_vector(fitted_values_unlisted)
  expect_length(fitted_values_unlisted, nrow(d))
})

test_that("calibration of metapred one-stage fixed effect models is ok", { 
  skip_on_cran()
  f <- y ~ x
  mp_fe <- metapred(d, "k", formula = f, scope = f, family = binomial, estFUN = glm, two.stage = F, 
                    perfFUN = list("bin.cal.int", "cal.slope"),
                    genFUN = list("rema"),
                    gen.of.perf = "factorial")
  expect_gte(gen(mp_fe, 1), 0)
  expect_lte(gen(mp_fe, 1), .01)

  expect_gte(gen(mp_fe, 2), 1)
  expect_lte(gen(mp_fe, 2), 1.1)
})

test_that("metapred can estimate one-stage random effects models", { 
  skip_on_cran()
  
  f_ri <- y ~ x + (1|k)
  
  cal_slope_glmer <- function(p, y, family = binomial, ...) {
    metamisc:::cal_slope(p, y, estFUN = "glm", family = family)
  }
  
  mp_ri <- metapred(d, "k", formula = f_ri, scope = f_ri, family = binomial, estFUN = lme4::glmer, two.stage = F, 
                    perfFUN = list("bin.cal.int", cal_slope_glmer),
                    genFUN = list("rema"),
                    gen.of.perf = "factorial")
  expect_is(mp_ri, "metapred")
  
  expect_gte(gen(mp_ri, 1), 0)
  expect_lte(gen(mp_ri, 1), .01)
  
  expect_gte(gen(mp_ri, 2), 1)
  expect_lte(gen(mp_ri, 2), 1.1)
  
  f_re <- y ~ x + (1 + x | k)
  mp_re <- metapred(d, "k", formula = f_re, scope = f_re, family = binomial, estFUN = lme4::glmer, two.stage = F)
  expect_is(mp_re, "metapred") 
})
