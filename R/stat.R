### stat functions
# extractAIC.crr, AIC.crr, BIC.crr, extractIC, logLik.crr, deviance.crr,
# coef.crr, coefficients.crr, terms.crr2, terms.crr, crrFits, crrWald.test
# 
# unexported: wald.test
###


#' \code{crr} model fit statistics
#' 
#' @description
#' Functions to assess the quality of fitted \code{\link[cmprsk]{crr}}
#' and \code{\link{crr2}} objects.
#' 
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
#' @param object an object of class \code{\link[cmprsk]{crr}}
#' @param ic information criterion, one of \code{"AIC"}, \code{"BIC"},
#' \code{"AICc"}, or \code{"BICc"}; see details
#' @param p an optional penalty term to be multiplied by, \code{k}, the
#' number of free parameters estimated in each model including the intercept
#' term
#' @param fit,scale,k see \code{\link{extractAIC}}
#' @param ... additional arguments passed to or from other methods
#' 
#' @seealso
#' \code{\link{crrFits}}; \code{\link{crrwald.test}}
#' 
#' @references
#' Kuk D, Varadhan R. Model selection in competing risks regression.
#' \emph{Stat Med}. 2013 Aug \strong{15};32.
#' 
#' @examples
#' crrs <- crr2(Surv(futime, event(censored) == death) ~ age, transplant)
#' crr1 <- crrs[[1L]]
#' 
#' crr1$coef
#' coef(crr1)
#' coefficients(crr1)
#' 
#' extractAIC(crr1)
#' sapply(crrs, AIC)
#' 
#' ## these are equivalent
#' extractIC(crrs[[1L]], p = 0)
#' -2 * logLik(crrs[[1]])[1]
#' -2 * crrs[[1]]$loglik
#' 
#' deviance(crr1)
#' 
#' crrFits(crr1)
#' crrwald.test(crr1)
#' 
#' @name crrfit
NULL

#' @rdname crrfit
#' @export
extractIC <- function(object, ic = c('AIC', 'BIC', 'AICc', 'BICc'), p) {
  cl <- object$call
  ic <- match.arg(ic)
  n  <- object[['n']]
  k  <- length(object[['coef']])
  
  ## get number of events of interest from function call
  ne <- tryCatch(
    sum(object$nuftime %||%
          eval(cl$fstatus) %in% eval(cl$failcode)),
    error = function(e) {
      if (ic == 'BIC')
        warning('BIC is approximate: use p = # of events of interest',
                call. = FALSE)
      length(object[['uftime']])
    })
  
  p <- if (!missing(p)) {
    ic <- paste('p', p, sep = '=')
    p
  } else switch(ic, AIC = 2, BIC = log(ne), AICc = 2,
                BICc = log(length(object[['res']][, 1L])))
  
  res <- p * k - 2 * object[['loglik']]
  
  if (ic == 'AICc')
    res <- res + 2 * (k + 1) * (k + 2) / (n - k - 2)
  setNames(res, ic)
}

#' @rdname crrfit
#' @export
extractAIC.crr <- function(fit, scale, k = 2, ...) {
  extractIC(fit, p = k)
}

#' @rdname crrfit
#' @export
AIC.crr <- function(object, ...) {
  extractIC(object, 'AIC')
}

#' @rdname crrfit
#' @export
BIC.crr <- function(object, ...) {
  extractIC(object, 'BIC')
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
  cat('Call:\n')
  dput(object$call)
  cat('\n')
  x  <- summary(object)
  dv <- unname(x$logtest[['test']])
  df <- x$logtest[['df']]
  pv <- pchisq(dv, df, lower.tail = FALSE)
  cat('Deviance =', signif(dv, 3L), 'on', df, 'df,',
      format.pval(pv, show.p = TRUE), '\n\n')
  invisible(c(chi.sq = dv, df = df, p.value = pv))
}

#' \code{crr} model coefficients
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

#' \code{crr2} model terms
#' 
#' @param x an object of class \code{\link{crr2}}
#' @param ... ignored
#' 
#' @export

terms.crr2 <- function(x, ...) {
  if (!inherits(x, 'crr2')) {
    message('\'terms\' is not available for a \'crr\' object - try \'?crr2\'')
    return(invisible(NULL))
  }
  terms(attr(x, 'model.frame'))
}

#' @rdname terms.crr2
#' @export
terms.crr <- terms.crr2

#' \code{crr} model selection table
#' 
#' Return several types of fit statistics for a \code{\link[cmprsk]{crr}}
#' saturated model compared to the null model. Note that comparisons among
#' models only make sense for ones fit to the same data set.
#' 
#' @param ... one or more objects of class \code{\link[cmprsk]{crr}}
#' @param p an optional penalty term to be multiplied by, \code{k}, the
#' number of free parameters estimated in each model including the intercept
#' term
#' 
#' @seealso
#' \code{\link{crrfit}}; \code{\link{crrwald.test}}
#' 
#' @references
#' Scrucca L, Santucci A, Aversa F (2009). Regression Modeling of Competing
#' Risk Using R: An In Depth Guide for Clinicians. \emph{Bone Marrow
#' Transplantation} (2010) \strong{45}, 1388-1395.
#' 
#' @examples
#' \dontrun{
#' ## example, figures, tables from
#' ## http://www.nature.com/bmt/journal/v45/n9/full/bmt2009359a.html
#' 
#' bmt <- read.csv('http://www.stat.unipg.it/luca/R/bmtcrr.csv')
#' bmt <- within(bmt, {
#'   Sex <- relevel(Sex, 'M')
#'   Phase <- relevel(Phase, 'Relapse')
#' })
#' cov1 <- with(bmt, model.matrix(~ Age + Sex + D + Phase + Source)[, -1])
#' 
#' m1 <- with(bmt, crr(ftime, Status, cov1))
#' summary(m1)
#' crrwald.test(m1, c(Phase = 'Phase'))
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
#' par(mfrow = c(2,2))
#' with(m2, {
#'   for (ii in 1:ncol(res))
#'     scatter.smooth(uftime, res[, ii],
#'                    main = gsub('Phase', '', names(coef)[ii]),
#'                    xlab = 'Failure time',
#'                    ylab = 'Schoenfeld residuals')
#' })
#' 
#' pred <- with(bmt, model.matrix(~ levels(Phase)))[, -1L]
#' pred <- predict(m2, pred)
#' 
#' plot(pred, xlab = 'Failure time', ylab = 'CIF', col = 1:4, lty = 1,
#'      ylim = c(0, 1))
#' legend('top', lty = 1, col = 1:4, horiz = TRUE, bty = 'n',
#'        title = 'Phase', legend = levels(bmt$Phase),
#'        x.intersp = .1, y.intersp = 0.5)
#' }
#' 
#' @export

crrFits <- function(..., p) {
  l <- list(...)
  
  ## crr2 objects are fit to same data and inherit crr already
  if (inherits(..1, 'crr2')) {
    l <- if (is.null(coef(..1)))
      c(...) else list(..1)
    ## drop cox models
    l <- Filter(function(x) !inherits(x, 'coxph'), l)
  } else
    ## for a manual list of crr objects?
    stopifnot(
      all(sapply(l, inherits, 'crr')),
      ## assert same data was used to fit
      # all(diff(sapply(l, function(x) x[['loglik.null']])) == 0),
      all(diff(sapply(l, function(x) x[['n']])) == 0)
    )
  
  null <- l[[1L]]
  null$loglik <- null$loglik.null
  null$coef <- null$call$cov1 <- null$call$cov2 <- NULL
  
  l <- c(list(null), l)
  n <- length(l)
  # model <- c('Null model', sapply(l[-1L], function(x) x[['call']]))
  
  model <- sapply(seq_along(l), function(ii)
    capture.output({
      if (ii == 1L)
        cat('0: Null Model', '\n')
      else {
        cat(sprintf('%s: Model %1$s call:\n', ii - 1L))
        dput(l[[ii]][['call']])
      }
      cat('\n')
    }))
  
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
    res[[sprintf('p=%s', p)]] <- pv <- vapply(l, extractIC, p = p, numeric(1L))
    res[[sprintf('p=%s diff', p)]] <- pv - min(pv)
  }
  
  res <- as.data.frame(res, check.names = FALSE)
  rownames(res) <- seq.int(n) - 1L
  
  title <- 'Model selection table\n'
  # model <- sprintf('Model %s: %s', rownames(res), model)
  
  structure(
    res, heading = c(title, unlist(model), ''),
    class = c('anova', 'data.frame')
  )
}

#' Wald test for \code{crr} model coefficients
#' 
#' @param object an object of class \code{\link[cmprsk]{crr}}
#' @param terms a list (optionally named) giving the indices of the
#' coefficients to test; a list of length \code{n} will perform \code{n}
#' tests; default is an overall test and each coefficient individually
#' 
#' Alternatively, a (named) vector (or list) of character strings can be
#' used which will \code{\link{grep}} for each pattern in the coefficient
#' names; patterns that match multiple terms will be tested as a group, and
#' only one character string per test will be used
#' 
#' @seealso
#' \code{\link[aod]{wald.test}}; \code{\link{crrfit}}; \code{\link{crrFits}}
#' 
#' @examples
#' crrs <- crr2(Surv(futime, event(censored) == death) ~
#'                age + sex + abo, transplant)
#' 
#' ## default is to test all estimates overall and individually
#' crrwald.test(crrs[[1L]])
#' summary(crrs[[1]])$coef[, 'p-value']
#' 
#' 
#' ## identical ways to call crrwald.test
#' crrwald.test(crrs[[1L]], list(overall = 1:5, Age = 1,
#'                               Sex = 2, ABO = 3:5))
#' 
#' crrwald.test(crrs[[1L]], c(overall = '.', Age = 'age', ## or list()
#'                            Sex = 'sex', ABO = 'abo'))
#' 
#' @export

crrwald.test <- function(object, terms) {
  assert_class(object, 'crr')
  co <- coef(object)
  nn <- names(co)
  nt <- seq_along(co)
  
  if (missing(terms)) {
    terms <- c(list(nt), nt)
    names(terms) <- c('Overall', nn)
  }
  if (is.character(unlist(terms)))
    terms <- lapply(terms, function(x) grep(x, nn))
  
  wald <- lapply(terms, function(x)
    wald.test(object[['var']], co, x)$result$chi2)
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
  
  structure(
    list(Sigma = Sigma, b = b, Terms = Terms, H0 = H0,
         L = L, result = res, verbose = verbose, df = df),
    class = 'wald.test'
  )
}
