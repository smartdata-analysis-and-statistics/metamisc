context("metapred 6. One-stage.")

set.seed(1234)
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



test_that("metapred can estimate one-stage random effects models", { 
  skip_on_cran()
  
  f <- y ~ x + (1|k)
  
  mp_ri <- metapred(d, "k", formula = f, scope = f, family = binomial, estFUN = lme4::glmer, two.stage = F)
  expect_is(mp_ri, "metapred")
  
  f <- y ~ x + (1 + x | k)
  mp_re <- metapred(d, "k", formula = f, scope = f, family = binomial, estFUN = lme4::glmer, two.stage = F)
  expect_is(mp_re, "metapred") 
})
