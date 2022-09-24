# formula to term labels
# utility function for obtaining terms (predictors, or interactions etc) of formula
# f formula
# returns term.labels (predictor names)
f2tl <- function(f) 
  attr(terms(f), "term.labels")

# Add or remove predictors
# f formula
# x predictors to be added or removed.
# NOTE: if one predictor is to be added, and one to be removed, an error is raised.
# Use two calls: updateFormula(updateFormula(f, xremove), xadd) or the other way around.
updateFormula <- function(f, x) {
  m <- x %in% f2tl(f)
  if (any(m)) if (!all(m)) stop("formula f must either contain all or no predictors x") # if any, all x must be contained in f 
  update.formula(f, paste("~ .", paste(if (any(m)) "-" else "+", as.character.default(x) ,
                                       sep = " ", collapse = " "), sep = " "  ) )
}

# Differences in vectors of characters
# x, y vector
# returns ordered vector of symmetric difference between x and y
setSymDiff <- function(x, y) 
  sort(unique(c(setdiff(x, y), setdiff(y, x))))

# Differences in formulas
# f, g formula
# returns ordered vector of symmetric difference between f and g 
getFormulaDiffAsChar <- function(f, g)
  setSymDiff(f2tl(f), f2tl(g))

# Get outcome of a formula
# f formula
# returns character, name of outcome
f2o <- function(f)
  as.character.default(f)[2]

# Get intercept-only formula from a formula
# f formula, containing outcome
# Returns formula, with only intercept on the right hand side.
f2iof <- function(f)
  formula(paste(f2o(f), "~ 1"))

# Get right hand side of formula from a formula
# f formula
# returns formula as: ~ right + hand + side
f2rhsf <- function(f)
  formula(paste(as.character.default(f)[1], as.character.default(f)[3]))

# Add two formulas up: sum of two formulas # NOTE: does not handle - 1
# f formula
# g formula
# Returns formula with left hand side of f, and both right hand sides.
addg2f <- function(f, g, terms = NULL) 
  formula(paste(f2o(f), "~", paste(c(union(f2tl(f), f2tl(g)), terms), collapse = " + ")))

# f <- y ~ x + z
# g <- y ~ x + a + I(a^2)
# h <- q ~ x + b - 1 + a*z
# 
# metamisc:::getFormulaDiffAsChar(f, g)
# metamisc:::f2o(f)
# metamisc:::addg2f(f, g)
# metamisc:::addg2f(f, h)
# metamisc:::addg2f(f, h, "studyid")