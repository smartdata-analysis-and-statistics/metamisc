## Note by Valentijn: I can't seem to make it pass all the checks for  forest.numeric or forest.default 
# (by Thomas, below) work without making a generic method. But that masks the generic from metafor. 
# The issue is that i was trying to unmask the generic from metafor...

#' Forest plot
#' 
#' Generate a forest plot by specifying the various effect sizes, confidence intervals and summary estimate.
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' @author Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>
#' 
#' @param \dots Additional arguments, which are currently ignored.
#' 
#' @details This is a generic function. See \link{forest.default} for making forest plots of summary statistics,
#' \link{forest.metapred} for plotting \link{metapred} objects, and \link{forest.mp.cv.val} for plotting 
#' \code{mp.cv.val} objects. 
#' 
#' @export forest
forest <- function(...)
  UseMethod("forest")

#' Forest plot
#' 
#' Generate a forest plot by specifying the various effect sizes, confidence intervals and summary estimate.
#' @param theta Numeric vector with effect size for each study
#' @param theta.ci.lb Numeric vector specifying the lower bound of the confidence interval of the effect sizes
#' @param theta.ci.ub Numeric vector specifying the upper bound of the confidence interval of the effect sizes
#' @param theta.slab Character vector specifying the study labels
#' @param theta.summary Meta-analysis summary estimate of the effect sizes
#' @param theta.summary.ci.lb Lower bound of the confidence (or credibility) interval of the summary estimate
#' @param theta.summary.ci.ub Upper bound of the confidence (or credibility) interval of the summary estimate
#' @param theta.summary.pi.lb Lower bound of the (approximate) prediction interval of the summary estimate. 
#' @param theta.summary.pi.ub Upper bound of the (approximate) prediction interval of the summary estimate.
#' @param title Title of the forest plot
#' @param sort By default, studies are sorted by ascending effect size (\code{sort="asc"}). Set to \code{"desc"} for 
#' sorting in reverse order, or any other value to ignore sorting.
#' @param theme Theme to generate the forest plot. By default, the classic dark-on-light ggplot2 theme is used. 
#' See \link[ggplot2]{ggtheme} for more information.
#' @param predint.linetype The linetype of the prediction interval
#' @param xlim The \code{x} limits \code{(x1, x2)} of the forest plot
#' @param xlab Optional character string specifying the X label
#' @param refline Optional numeric specifying a reference line
#' @param label.summary Optional character string specifying the label for the summary estimate
#' @param label.predint Optional character string specifying the label for the (approximate) prediction interval
#' @param nrows.before.summary How many empty rows should be introduced between the study results and the summary estimates
#' @param study.digits How many significant digits should be used to print the stuy results
#' @param study.shape Plotting symbol to use for the study results. By default, a filled square is used.
#' @param col.diamond The filling color for the diamond representing the summary estimate. 
#' E.g. "red", "blue", or hex color code ("#2e8aff")
#' @param col.predint Line color for the prediction interval. E.g. "red", "blue", or hex color code ("#2e8aff")
#' @param size.study Line width for the study results in mm
#' @param size.predint Line width for the prediction interval in mm
#' @param lty.ref Line type for the reference line
#' @param \dots Additional arguments, which are currently ignored.
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @method forest default
#' 
#' @return An object of class \code{ggplot}
#' 
#' @import metafor
#' @import ggplot2
#' @importFrom stats reorder
#' 
#' @export
# This could also be named forest.numeric
forest.default <- function(theta, 
                    theta.ci.lb,
                    theta.ci.ub,
                    theta.slab, 
                    theta.summary, 
                    theta.summary.ci.lb, theta.summary.ci.ub, 
                    theta.summary.pi.lb, theta.summary.pi.ub,
                    title,
                    sort = "asc",
                    theme = theme_bw(),
                    predint.linetype = 1,
                    xlim,
                    xlab = "", 
                    refline = 0,
                    label.summary = "Summary Estimate", 
                    label.predint = "Prediction Interval",
                    nrows.before.summary = 1,
                    study.digits = 2,
                    study.shape = 15,
                    col.diamond = "white",
                    col.predint = "black",
                    size.study = 0.5,
                    size.predint = 1,
                    lty.ref = "dotted",
                    ...) {
  requireNamespace("ggplot2")

  if (missing(theta)) stop("Study effect sizes are missing!")
  if (missing(theta.ci.lb) | missing(theta.ci.ub)) stop("Confidence intervals of effect sizes missing!")
  if (missing(theta.slab)) stop("Study labels are missing!")
  
  num.studies <- unique(c(length(theta), length(theta.ci.lb), length(theta.ci.ub), length(theta.slab)))
  if (length(num.studies) > 1) stop(paste("Too few studies for a forest plot!"))

  #Extract data
  yi <- theta
  k <- length(theta)
  
  if (missing(theta.slab)) {
    if (!is.null(attr(theta, "slab"))) {
      slab <- attr(theta, "slab")
    } else {
      slab <- paste("Study", seq_len(k))
    }
  } else {
    if (length(theta.slab) == 1 && is.na(theta.slab)) 
      slab <- rep("", k)
    else
      slab <- theta.slab
  }
  
  add.predint <- TRUE # Add prediction interval by default
  if (missing(theta.summary.pi.lb) | missing(theta.summary.pi.ub)) {
    theta.summary.pi.lb <- theta.summary.pi.ub <- NA
  }
  if (NA %in% c(theta.summary.pi.lb, theta.summary.pi.ub)) {
    add.predint <- FALSE
  }
  
  # Determine ordering
  i.index <- 1:length(yi)
  if (sort == "asc") {
    i.index <- order(yi, decreasing = FALSE)
  } else if (sort == "desc") {
    i.index <- order(yi, decreasing = TRUE)
  }
  
  # Order the study results
  scat  <- rep(1, num.studies) #indicator variable for data points
  slab  <- slab[i.index]
  yi    <- yi[i.index]
  ci.lb <- theta.ci.lb[i.index]
  ci.ub <- theta.ci.ub[i.index]
  
  if (!missing(theta.summary)) {
    if (!add.predint) {
      scat  <- c(scat, 0)
      slab  <- c(slab, label.summary)
      yi    <- c(yi, theta.summary)
      ci.lb <- c(ci.lb, theta.summary.ci.lb)
      ci.ub <- c(ci.ub, theta.summary.ci.ub)
    } else {
      scat  <- c(scat, 0, 0)
      slab  <- c(slab, label.summary, label.predint)
      yi    <- c(yi, theta.summary, theta.summary)
      ci.lb <- c(ci.lb, theta.summary.ci.lb, theta.summary.pi.lb)
      ci.ub <- c(ci.ub, theta.summary.ci.ub, theta.summary.pi.ub)
    }
  } 
  
  
  ALL <- data.frame("study" = slab, 
                    "mean" = yi, 
                    "m.lower" = ci.lb, 
                    "m.upper" = ci.ub, 
                    "order" = length(yi):1, 
                    "scat" = scat)
  
  # Add extra space between the study results and the summary estimates
  rows.summaries <- which(ALL$study %in% c(label.summary, label.predint))
  ALL$order[rows.summaries] <- ALL$order[rows.summaries] - nrows.before.summary

  # reorder factor levels based on another variable (by yi)
  ALL$study.ES_order <- reorder(ALL$study, ALL$order, mean) 
  
  # Information for secondary axis
  labels_axis2 <- paste(format(ALL$mean, digits = study.digits, nsmall=2), " [", 
                        format(ALL$m.lower, digits = study.digits, nsmall=2), " ; ", 
                        format(ALL$m.upper, digits = study.digits, nsmall=2), "]", sep = "")
  
  
  p <- with(ALL, ggplot(ALL[!is.na(ALL$mean), ], 
                        aes(x = order, y = mean, ymin = m.lower, ymax = m.upper)) +
              theme +
              theme(axis.ticks.y.right = element_blank(),
                    panel.grid.minor.y = element_blank()) +
              geom_pointrange(data = subset(ALL, scat == 1), shape = study.shape, size = size.study) + 
              scale_x_continuous(
                breaks   = order,
                labels   = study,
                sec.axis = ggplot2::sec_axis(
                  ~ .,
                  breaks = order,
                  labels = labels_axis2)
                ) +
              coord_flip() + 
              ylab(xlab) + 
              xlab(""))

  
  if (!missing(xlim)) {
    p <- p + ylim(xlim)
  }
  
  # Add title
  if (!missing(title)) {
    p <- p + ggtitle(title)
  }
  
  # Add refline
  if (!missing(refline)) {
    if (is.numeric(refline)) {
      p <- p + geom_hline(yintercept = refline,  linetype = lty.ref) 
    }
  }
  
  # Add meta-analysis summary
  if (!missing(theta.summary)) {
    g2 <- with(ALL, subset(ALL, study == label.summary))
    g2$ci.upper <- theta.summary.ci.ub
    g2$ci.lower <- theta.summary.ci.lb
    
    
    # Add confidence interval of the summary estimate
    p <- p + with(g2, geom_errorbar(data = g2, aes(ymin = ci.lower, ymax = ci.upper, x = order), width = 0.5, size=1.0))
    
    # Add summary estimate
    p <- p + with(g2, geom_point(data = g2, aes(x = order, y = mean), shape=23, size=3, fill = col.diamond))
    
    # Add diaomond
    #df_diamond <- data.frame(row = c(g2$order+0.5, g2$order, g2$order-0.5, g2$order), ci = c(0, 1, 0, 0))
    #print(df_diamond)
    #p <- p + geom_polygon(data = df_diamond, aes(x = row, y = ci))
    
    # Add (approximate) prediction interval
    if (add.predint) {
      g3 <- with(ALL, subset(ALL, study == label.predint))
      g3$pi.upper <- theta.summary.pi.ub
      g3$pi.lower <- theta.summary.pi.lb
      
      p <- p + with(g3, geom_errorbar(data = g3, 
                             aes(ymin = pi.lower, ymax = pi.upper, x = order), 
                             size = size.predint, 
                             width = 0.5,
                             color = col.predint,
                             linetype = predint.linetype))
    }
  }
  

  p
}

#' Forest Plots
#' 
#' Function to create forest plots for objects of class \code{"mm_perf"}.
#' 
#' @param x An object of class \code{"mm_perf"}
#' @param \ldots Additional arguments which are passed to \link{forest}.
#' 
#' @details The forest plot shows the performance estimates of each study with corresponding confidence 
#' intervals. 
#' 
#' @references 
#' Lewis S, Clarke M. Forest plots: trying to see the wood and the trees. \emph{BMJ}. 2001; 322(7300):1479--80.
#' 
#' @examples 
#' data(EuroSCORE)
#' 
#' # Calculate the c-statistic and its standard error
#' est1 <- ccalc(cstat = c.index, cstat.se = se.c.index, cstat.cilb = c.index.95CIl, 
#'               cstat.ciub = c.index.95CIu, N = n, O = n.events, data = EuroSCORE, slab = Study)
#' plot(est1)
#' 
#' # Calculate the total O:E ratio and its standard error
#' est2 <- oecalc(O = n.events, E = e.events, N = n, data = EuroSCORE, slab = Study)
#' plot(est2)
#' 
#' @keywords forest
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @method plot mm_perf
#' @export
plot.mm_perf <- function(x, ...) {
  forest(theta = x$theta, theta.ci.lb = x$theta.cilb, theta.ci.ub = x$theta.ciub, theta.slab = rownames(x), ...)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  #library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Posterior distribution of estimated model parameters
#' 
#' Generate a plot of the posterior distribution
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot} object.
#' 
#' @details This is a generic function. 
#' 
#' @export dplot
dplot <- function(...)
  UseMethod("dplot")

#' Gelman-Rubin-Brooks plot
#' 
#' This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases. The code is adapted from
#' the R package coda.
#' 
#' @param \ldots Additional arguments which are currently not used
#' @return A \code{ggplot} object.
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @export gelmanplot
gelmanplot <- function(...) 
  UseMethod("gelmanplot")

#' Plot the autocorrelation of a Bayesian meta-analysis model
#' 
#' Function to display autocorrelation of a fitted Bayesian meta-analysis model.
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot} object.
#'
#' @details This is a generic function. 
#' 
#' @export acplot
acplot <- function(...)
  UseMethod("acplot")

#' Plot the running means of a Bayesian meta-analysis model
#' 
#' Function to display running means of a fitted Bayesian meta-analysis model.
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' @param \dots Additional arguments, which are currently ignored.
#' @return A \code{ggplot} object.
#'  
#' @details This is a generic function. 
#' 
#' @export rmplot
rmplot <- function(...)
  UseMethod("rmplot")

#' Plot the autocorrelation of a Bayesian meta-analysis model
#' 
#' Function to display autocorrelation of a fitted Bayesian meta-analysis model.
#' 
#' @param x An object of class \code{"mcmc.list"}
#' @param P Optional dataframe describing the parameters to plot and their respective names
#' @param greek Logical value indicating whether parameter labels have to be parsed to get Greek letters. Defaults to FALSE.
#' @param nLags	Integer indicating the number of lags of the autocorrelation plot.
#' @param \ldots Additional arguments which are passed to ggs_autocorrelation
#' @return A \code{ggplot} object.
#' 
#' @keywords meta-analysis autocorrelation
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
#' @importFrom dplyr group_by do %>%
acplot.mcmc.list <- function(x, P, nLags = 50, greek = FALSE, ...) {
  requireNamespace("ggmcmc")
  
  if (!missing(P)) {
    S <- ggmcmc::ggs(x, par_labels = P, sort = FALSE)
    D <- subset(S, S$ParameterOriginal %in% P$Parameter)
  } else {
    D <- ggmcmc::ggs(x, sort = FALSE)
  }
  g <- ggmcmc::ggs_autocorrelation(D = D, nLags = nLags, greek = greek, ...)
  g <- g + theme(strip.text.y = element_text(angle = 0))
  
  return(g)
}


#' Plot the running means of a Bayesian meta-analysis model
#' 
#' Function to display running means of a fitted Bayesian meta-analysis model.
#' 
#' @param x An object of class \code{"mcmc.list"}
#' @param P Optional dataframe describing the parameters to plot and their respective names
#' @param greek Logical value indicating whether parameter labels have to be parsed to get Greek letters. Defaults to FALSE.
#' @param \ldots Additional arguments which are passed to ggs_running
#' @return A \code{ggplot} object.
#' 
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
rmplot.mcmc.list <- function(x, P, greek = FALSE, ...) {
  requireNamespace("ggmcmc")
  
  if (!missing(P)) {
    S <- ggmcmc::ggs(x, par_labels = P, sort = FALSE)
    D <- subset(S, S$ParameterOriginal %in% P$Parameter)
  } else {
    D <- ggmcmc::ggs(x, sort = FALSE)
  }
  
  g <- ggmcmc::ggs_running(D = D, greek = greek, ...)
  g <- g + theme(strip.text.y = element_text(angle = 0))
  
  return(g)
}


#' Posterior distribution of estimated model parameters
#' 
#' Generate a plot of the posterior distribution
#' 
#' @param x An object of class \code{"mcmc.list"}
#' @param P Optional dataframe describing the parameters to plot and their respective names
#' @param plot_type Optional character string to specify whether a density plot (\code{"dens"}) or 
#' histogram (\code{"hist"}) should be displayed.
#' @param \ldots Additional arguments which are currently not used
#' @return A \code{ggplot} object.
#' 
#' 
#' @keywords meta-analysis density distribution
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
dplot.mcmc.list <- function(x, P, plot_type = "dens", ...) {
  requireNamespace("ggmcmc")
  
  if (!missing(P)) {
    S <- ggmcmc::ggs(x, par_labels = P, sort = FALSE)
    S <- subset(S, S$ParameterOriginal %in% P$Parameter)
  } else {
    S <- ggmcmc::ggs(x, sort = FALSE)
  }
  
  if (plot_type == "dens") {
    ggmcmc::ggs_density(S)
  } else if (plot_type == "hist") {
    ggmcmc::ggs_histogram(S)
  } else {
    stop("Invalid plot type")
  }
}

#' Gelman-Rubin-Brooks plot
#' 
#' This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases. The code is adapted from
#' the R package coda.
#' 
#' @param x An mcmc object
#' @param P Optional dataframe describing the parameters to plot and their respective names
#' @param confidence The coverage probability of the confidence interval for the potential scale reduction factor
#' @param max.bins Maximum number of bins, excluding the last one.
#' @param autoburnin Logical flag indicating whether only the second half of the series should be used in the computation. 
#' If set to TRUE (default) and start(x) is less than end(x)/2 then start of series will be adjusted so that only second half of series is used.
#' @param greek Logical value indicating whether parameter labels have to be parsed to get Greek letters. Defaults to false.
#' @param \ldots Additional arguments which are currently not used
#' @return A \code{ggplot} object.
#' 
#'             
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @return An object of class \code{ggplot}
#' 
#' @export
gelmanplot.mcmc.list <- function(x, P, confidence = 0.95, max.bins = 50, autoburnin = TRUE, greek = FALSE, ...) {
  requireNamespace("coda")

  nbin <- min(floor((coda::niter(x) - 50)/coda::thin(x)), max.bins)
  if (nbin < 1) {
    stop("Insufficient iterations to produce Gelman-Rubin plot")
  }
  binw <- floor((coda::niter(x) - 50)/nbin)
  last.iter <- c(seq(from = start(x) + 50 * coda::thin(x), by = binw * 
                       coda::thin(x), length = nbin), end(x))
  shrink <- array(dim = c(nbin + 1, coda::nvar(x), 2))
  dimnames(shrink) <- list(last.iter, coda::varnames(x), 
                           c("median", paste(50 * (confidence + 1), "%", sep = "")))
  for (i in 1:(nbin + 1)) {
    shrink[i, , ] <- coda::gelman.diag(window(x, end = last.iter[i]), 
                                 confidence = confidence, 
                                 transform = FALSE, 
                                 autoburnin = autoburnin, 
                                 multivariate = FALSE)$psrf
  }
  all.na <- apply(is.na(shrink[, , 1, drop = FALSE]), 2, all)
  if (any(all.na)) {
    cat("\n******* Error: *******\n")
    cat("Cannot compute Gelman & Rubin's diagnostic for any chain \n")
    cat("segments for variables", coda::varnames(x)[all.na], "\n")
    cat("This indicates convergence failure\n")
  }
  
  labels_gelmandiag <- names(shrink[1,1,])
  ggdatMed <- data.frame(cbind(shrinkfactor = as.vector(shrink[,,labels_gelmandiag[1]]), 
                               parameter = rep(colnames(shrink[,,1]), each = nrow(shrink[,,1])),
                               type = labels_gelmandiag[1], 
                               last.iter = rep(last.iter, ncol(shrink[,,1]))))
  ggdatBnd <- data.frame(cbind(shrinkfactor = as.vector(shrink[,,labels_gelmandiag[2]]), 
                               parameter = rep(colnames(shrink[,,2]), each = nrow(shrink[,,2])),
                               type = labels_gelmandiag[2], 
                               last.iter = rep(last.iter, ncol(shrink[,,2]))))
  ggdat <- rbind(ggdatMed, ggdatBnd)
  ggdat$last.iter <- as.numeric(as.character(ggdat$last.iter))
  ggdat$shrinkfactor <- as.numeric(as.character(ggdat$shrinkfactor))
  
  if (!missing(P)) {
    ggdat <- subset(ggdat, ggdat$parameter %in% P$Parameter)
    ggdat$parameter <- factor(ggdat$parameter, levels = P$Parameter, labels = P$Label)
  }
  
  
  if (greek) {
    g <- with(ggdat, ggplot(ggdat, aes(last.iter, shrinkfactor, color = type)) +
                geom_line() + facet_wrap(~ parameter, labeller = label_parsed) + 
                xlab("Last iteration in chain") +
                ylab("Shrink factor"))
  } else {
    g <- with(ggdat, ggplot(ggdat, aes(last.iter, shrinkfactor, color = type)) +
                geom_line() + facet_wrap(~ parameter) + 
                xlab("Last iteration in chain") +
                ylab("Shrink factor"))
  }
  
  
  return(g)
}