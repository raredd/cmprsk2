### stat
# crrfit, extractAIC.crr, AIC.crr, extractBIC.crr, BIC.crr, extractIC,
# logLik.crr, deviance.crr, coef.crr, coefficients.crr
###


### TODO
# terms.crr
###


#' \code{crr} model fit statistics
#' 
#' Functions to assess the quality of fitted \code{\link[cmprsk]{crr}} models.
#' 
#' @param object an object of class \code{\link[cmprsk]{crr}}
#' @param ... ignored
#' 
#' @references
#' Kuk D, Varadhan R. Model selection in competing risks regression.
#' \emph{Stat Med}. 2013 Aug \strong{15};32.
#' 
#' @name crrfit
NULL

#' Extract AIC from a \code{crr} model
#' 
#' @param k numeric specifying the "weight" of the equivalent degrees of
#' freedom (\code{=: edf}) part in the AIC formula
#' 
#' @examples
#' fit <- crr2(Surv(futime, event(censored) == death) ~ age, transplant)
#' extractAIC(fit[[1]])
#' AIC(fit[[1]])
#' 
#' @rdname crrfit
#' @export
extractAIC.crr <- function(object, ...) {
  extractIC(object, 'AIC')
}

#' @rdname crrfit
#' @export
AIC.crr <- extractAIC.crr

#' @rdname crrfit
#' @export
extractBIC.crr <- function(object, ...) {
  extractIC(object, 'BIC')
}

#' @rdname crrfit
#' @export
BIC.crr <- extractBIC.crr

#' @description
#' Select multiple types of Akaike or Bayesian (Schwarz) information
#' criterion from a \code{\link[cmprsk]{crr}} object to assess the relative
#' quality of models for a given data set.
#' 
#' @details
#' \code{AIC} and \code{BIC} are calculated in the usual way. \code{AICc}
#' is the AIC with a correction for finite sample sizes. This assumes that
#' the model is univariate, linear, and has normally-distributed residuals
#' (conditional upon regressors).
#' 
#' \code{BICc} is a newly proposed criteria that is a modified BIC for
#' competing risks data subject to right censoring (Kuk, 2013).
#' 
#' @param ic information criterion, one of \code{"AIC"}, \code{"BIC"},
#' \code{"AICc"}, or \code{"BICc"}; see details
#' 
#' @references
#' Kuk D, Varadhan R. Model selection in competing risks regression.
#' \emph{Stat Med}. 2013 Aug \strong{15};32.
#' 
#' @examples
#' 
#' @rdname crrfit
#' @export

extractIC <- function(object, ic = c('AIC', 'BIC', 'AICc', 'BICc'), p) {
  n  <- object[['n']]
  cl <- object$call
  
  ## get number of events of interest from function call
  ne <- tryCatch(
    sum(eval(cl$fstatus) %in% eval(cl$failcode)),
    error = function(e) {
      warning('BIC is approximate: use p = #events of interest', call. = FALSE)
      length(object[['uftime']])
    })
  
  k  <- length(object[['coef']])
  ic <- match.arg(ic)
  
  p <- if (!missing(p)) {
    ic <- 'p'
    p
  } else switch(ic, AIC = 2, BIC = log(ne), AICc = 2,
                BICc = log(length(object[['res']][, 1L])))
  
  res <- p * k - 2 * object[['loglik']]
  
  if (ic == 'AICc')
    res <- res + 2 * (k + 1) * (k + 2) / (n - k - 2)
  res
}

#' @rdname crrfit
#' @export
logLik.crr <- function(object, ...) {
  val <- object[['loglik']]
  attr(val, 'df') <- length(object[['coef']])
  class(val) <- 'logLik'
  val
}

#' @rdname crrfit
#' @export
deviance.crr <- function(object, ...) {
  x  <- summary(object)
  dv <- unname(x$logtest[['test']])
  df <- x$logtest[['df']]
  pv <- pchisq(dv, df, lower.tail = FALSE)
  cat('Deviance =', signif(dv, 3L), 'on', df, 'df,',
      format.pval(pv, show.p = TRUE), '\n\n')
  invisible(c(chi.sq = dv, df = df, p.value = pv))
}

#' ct \code{crr} model coefficients
#' 
#' @param object an object of class \code{\link[cmprsk]{crr}}
#' @param ... ignored
#' 
#' @export

coef.crr <- function(object, ...) {
  object[['coef']]
}

#' @rdname coef.crr
#' @export
coefficients.crr <- coef.crr

#' Model terms
#' 
#' @param x an object of class \code{\link[cmprsk]{crr}}
#' 
#' @export

terms.crr <- function(x, ...) {
  
}

#' \code{crr} model selection table
#' 
#' Return several types of fit statistics for a \code{\link[cmprsk]{crr}}
#' saturated model compared to the null model. Note that comparisons among
#' models only make sense for ones fit to the same data set.
#' 
#' @param ... one or more objects of class \code{\link[cmprsk]{crr}}
#' @param p a penalty term to be multiplied by, \code{k}, the number of free
#' parameters estimated in each model including the intercept term
#' 
#' @references
#' Scrucca L, Santucci A, Aversa F (2009). Regression Modeling of Competing
#' Risk Using R: An In Depth Guide for Clinicians. \emph{Bone Marrow
#' Transplantation} (2010) \strong{45}, 1388–1395.
#' 
#' @examples
#' bmt <- read.csv('http://www.stat.unipg.it/luca/R/bmtcrr.csv')
#' cov1 <- with(bmt, model.matrix(~ Age + relevel(Sex, 'M') + D +
#'              relevel(Phase, 'Relapse') + Source)[, -1])
#' m1 <- with(bmt, crr(ftime, Status, cov1))
#' summary(m1)
#' aod::wald.test(m1$var, m1$coef, grep('Phase', names(m1$coef)))$result$chi2
#' 
#' ## model selection
#' m2 <- with(bmt, crr(ftime, Status, cov1[, c(4:6)]))
#' m3 <- with(bmt, crr(ftime, Status, cov1[, c(4:6,7)]))
#' m4 <- with(bmt, crr(ftime, Status, cov1[, c(4:6,7,1)]))
#' m5 <- with(bmt, crr(ftime, Status, cov1[, c(4:6,7,2)]))
#' m6 <- with(bmt, crr(ftime, Status, cov1[, c(4:6,7,3)]))
#' 
#' crrFits(m1, m2, m3, m4, m5, m6, p = 3)
#' 
#' par(mfrow = c(1,3))
#' for (ii in 1:ncol(m2$res))
#'   scatter.smooth(m2$uftime, m2$res[, ii],
#'                  main = names(coef(m2))[ii],
#'                  xlab = 'Failure time',
#'                  ylab = 'Schoenfeld residuals')
#' 
#' pred <- with(bmt, model.matrix(~ levels(relevel(Phase, 'Relapse')))[, -1])
#' pred <- predict(m2, pred)
#' plot(pred, xlab = 'Failure time', ylab = 'CIF', lty = 1:4)
#' legend('topleft', lty = 1:4, title = 'Phase', bty = 'n',
#'        legend = levels(relevel(bmt$Phase, "Relapse")))
#' 
#' 
#' @export

crrFits <- function(..., p) {
  l <- list(...)
  stopifnot(
    all(sapply(l, inherits, 'crr')),
    ## assert same data was used to fit
    # all(diff(sapply(l, function(x) x[['loglik.null']])) == 0),
    all(diff(sapply(l, function(x) x[['n']])) == 0)
  )
  
  null <- ..1
  null$loglik <- null$loglik.null
  null$coef <- null$call$cov1 <- null$call$cov2 <- NULL
  
  l <- c(list(null), l)
  n <- length(l)
  model <- c('Null model', sapply(l[-1L], function(x) x[['call']]))
  
  res <- list()
  res$n      <- vapply(l, function(x) x$n, integer(1L))
  res$loglik <- vapply(l, logLik, numeric(1L))
  res$df     <- vapply(l, function(x) length(x$coef), integer(1L))
  res$k      <- res$df + 1L
  
  res$`-2logLik`      <- vapply(l, extractIC, p = 0, numeric(1L))
  res$`-2logLik diff` <- res$`-2logLik` - min(res$`-2logLik`)
  
  res$AIC        <- vapply(l, AIC, numeric(1L))
  res$`AIC diff` <- res$AIC - min(res$AIC)
  
  res$BIC        <- vapply(l, BIC, numeric(1L))
  res$`BIC diff` <- res$BIC - min(res$BIC)
  
  if (!missing(p)) {
    res$Criterion        <- vapply(l, extractIC, p = p, numeric(1L))
    res$`Criterion diff` <- res$Criterion - min(res$Criterion)
  }
  
  res <- as.data.frame(res, check.names = FALSE)
  rownames(res) <- seq.int(n) - 1L
  
  title <- 'Model selection table\n'
  model <- sprintf('Model %s: %s', rownames(res), model)
  
  structure(res, heading = c(title, model, ''),
            class = c('anova', 'data.frame'))
}

#' Wald test for \code{crr} model coefficients
#' 
#' @param object an object of class \code{\link[cmprsk]{crr}}
#' @param terms a list (optionally named) giving the coefficients to test;
#' default is an overall test and each coefficient individually
#' 
#' @seealso
#' \code{\link[aod]{wald.test}}
#' 
#' @examples
#' fits <- crr2(Surv(t2, status(0) == 1) ~ Group + logtime, bmt2)
#' crrWald(fits[[1L]], list(overall = 1:3, Group = 1:2, logtime = 3))
#' 
#' @export

crrWald <- function(object, terms) {
  stopifnot(inherits(object, 'crr'))
  nt <- seq.int(length(object[['coef']]))
  if (missing(terms))
    terms <- c(list(nt), nt)
  wald <- lapply(terms, function(x)
    wald.test(object[['var']], object[['coef']], x)$result$chi2)
  res <- do.call('rbind', wald)
  rownames(res) <- names(terms)
  res
}

## aod::wald.test
wald.test <- function (Sigma, b, Terms = NULL, L = NULL, H0 = NULL,
                       df = NULL, verbose = FALSE) {
  if (is.null(Terms) & is.null(L)) 
    stop("One of the arguments Terms or L must be used.")
  if (!is.null(Terms) & !is.null(L)) 
    stop("Only one of the arguments Terms or L must be used.")
  if (is.null(Terms)) {
    w <- nrow(L)
    Terms <- seq(length(b))[colSums(L) > 0]
  } else w <- length(Terms)
  if (is.null(H0)) 
    H0 <- rep(0, w)
  if (w != length(H0)) 
    stop("Vectors of tested coefficients and of null hypothesis have different lengths\n")
  if (is.null(L)) {
    L <- matrix(rep(0, length(b) * w), ncol = length(b))
    for (i in 1:w)
      L[i, Terms[i]] <- 1
  }
  dimnames(L) <- list(paste0('L', as.character(seq(NROW(L)))), names(b))
  f <- L %*% b
  V <- Sigma
  mat <- qr.solve(L %*% V %*% t(L))
  stat <- t(f - H0) %*% mat %*% (f - H0)
  p <- 1 - pchisq(stat, df = w)
  if (is.null(df)) 
    res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  else {
    fstat <- stat/nrow(L)
    df1 <- nrow(L)
    df2 <- df
    res <- list(chi2 = c(chi2 = stat, df = w, P = p),
                Ftest = c(Fstat = fstat, df1 = df1, df2 = df2,
                          P = 1 - pf(fstat, df1, df2)))
  }
  structure(list(Sigma = Sigma, b = b, Terms = Terms, H0 = H0, 
                 L = L, result = res, verbose = verbose, df = df),
            class = 'wald.test')
}
