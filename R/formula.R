### utils for working with formulas
# parse_formula, sterms.inner, strata, terms.inner
# 
# formula s3 methods:
# is.crr2, is.crr2.default, is.crr2.formula, is.cuminc2, is.cuminc2.default,
# is.cuminc2.formula
###


is.crr2 <- function(x) {
  ## method to validate crr2 objects, formulae
  if (inherits(x, 'crr2'))
    return(TRUE)
  UseMethod('is.crr2')
}

is.crr2.default <- function(x) {
  inherits(x, 'crr2')
}

is.crr2.formula <- function(x) {
  x <- if (is.character(x) ||
           length(rapply(as.list(x)[2L], function(y)
             length(as.list(y)))) == 1L)
    deparse1(x) else deparse1(x[[2L]])
  
  ## checks for status(0) but not status(0) == 1
  # grepl('[^(]+\\([^(~]+\\((?=.+~|[^~]+$)', x, perl = TRUE)
  
  ## asserts status(0) == 1
  grepl('Surv\\([^(]+\\([^(]+\\)\\s*(==|%in%)[^)]+\\)', x)
}

is.cuminc2 <- function(x) {
  ## method to validate cuminc2 objects, formulae
  if (inherits(x, 'cuminc2'))
    return(TRUE)
  UseMethod('is.cuminc2')
}

is.cuminc2.default <- function(x) {
  inherits(x, 'cuminc2')
}

is.cuminc2.formula <- function(x) {
  x <- if (is.character(x) ||
           length(rapply(as.list(x)[2L], function(y)
             length(as.list(y)))) == 1L)
    deparse1(x) else deparse1(x[[2L]])
  
  ## status(0) == 1 equality is optional
  grepl('Surv\\([^(]+\\([^(]+\\)\\s*(?:(==|%in%)[^)]+)?\\)', x) &
    grepl('~', x)
}

#' Parse \code{cmprsk2} formula
#' 
#' Parse a \code{\link{crr2}} or \code{\link{cuminc2}} formula into useful
#' pieces and optionally check data source for errors.
#' 
#' @param formula a valid \code{crr2} or \code{cuminc2} formula
#' @param data,name optional data frame and data name to check variables
#'   from \code{formula}
#' 
#' @return
#' A list with the following elements:
#' 
#' \item{\code{formula}}{the unchanged \code{formula} coerced to a
#'   \code{\link{formula}} object}
#' \item{\code{lhs_vars}}{all variable names found on the left-hand side of
#'   \code{formula}}
#' \item{\code{rhs_vars}}{all variable names found on the right-hand side of
#'   \code{formula} with any \code{factor()}, \code{strata()}, etc.
#'   modifications removed}
#' \item{\code{lhs}}{the left-hand side of \code{formula} as it appears}
#' \item{\code{rhs}}{the right-hand side of \code{formula} as it appears with
#'   any \code{factor()}, \code{strata()}, etc. included}
#' \item{\code{ftime}}{variable in \code{data} defining times}
#' \item{\code{fstatus}}{variable in \code{data} defining censoring and
#'   failure events}
#' \item{\code{cencode}}{censor code, a value of \code{data$fstatus}}
#' \item{\code{failcode}}{failure code of interest if given, a value of
#'   \code{[data$fstatus]} (required for a \code{crr2} model but optional for
#'   \code{cuminc2})}
#' \item{\code{cov1}}{\code{cov1} variable names to be evaluated by
#'   \code{\link{model.matrix}} as they appear}
#' \item{\code{cov2}}{\code{cov2} variable names (functions of time) to be
#'   evaluated by \code{\link{model.matrix}} as they appear}
#' \item{\code{strata}}{strata variable names}

parse_formula <- function(formula, data = NULL, name = NULL) {
  ## parse crr2/cuminc2 formulae into useful pieces
  # cmprsk2:::parse_formula(Surv(time, status(censor) == event) ~ x)
  # cmprsk2:::parse_formula(Surv(time, status(censor) == event) ~ x + tf(y) + tf(factor(z)))
  formula <- as.formula(formula)
  
  if (!(is.crr2(formula) | is.cuminc2(formula)))
    stop('Invalid formula - see ?crr2 or ?cuminc2 for details')
  
  tt <- terms(formula, data = data, specials = c('tf', 'strata'))
  vv <- rownames(attr(tt, 'factors'))
  sp <- attr(tt, 'specials')
  ti <- terms.inner(formula)
  
  lhs_vars <- ti[[1L]][1:2]
  rhs_vars <- if (is.null(vv)) NULL else all.vars(reformulate(vv[-1L]))
  strata   <- strata(formula)
  ftime    <- ti[[1L]][1L]
  fstatus  <- ti[[1L]][2L]
  cencode  <- ti[[2L]][[1L]] ## required
  failcode <- sort(ti[[2L]][2L]) %||% NULL ## possible NA if not given, ie, cuminc2
  
  if (!is.null(data)) {
    dname <- name %||% deparse1(substitute(data))
    vars  <- c(lhs_vars, rhs_vars, strata)
    
    if (length(vars <- vars[!vars %in% names(data)]))
      stop(
        sprintf('%s not found in %s', toString(shQuote(vars)), shQuote(dname))
      )
    
    if (!is.numeric(na.omit(data[, ftime, drop = TRUE])) |
        any(na.omit(data[, ftime, drop = TRUE]) < 0))
      stop(
        sprintf('%s should be numeric values >= 0', shQuote(ftime))
      )
    
    codes <- na.omit(c(cencode, failcode))
    if (length(codes <- codes[!codes %in% data[, fstatus, drop = TRUE]]))
      warning(
        sprintf('%s not found in %s[, %s]',
                toString(shQuote(codes)), dname, shQuote(fstatus))
      )
    
    if (length(wh <- unique(data[, fstatus, drop = TRUE])) <= 2L)
      warning(
        gsub('\\n\\s{2,}', ' ',
             sprintf('Only %s found in in %s[, %s]\n Typically, censoring,
                   event of interest, and >=1 competing events are used',
                     toString(shQuote(wh)), deparse1(name), shQuote(fstatus)))
      )
  }
  
  list(
    formula  = formula,
    lhs_vars = lhs_vars,
    rhs_vars = rhs_vars,
    
    lhs = vv[1L],
    rhs = vv[-1L],
    
    ftime    = ftime,
    fstatus  = fstatus,
    cencode  = cencode,
    failcode = failcode,
    
    ## cov1 ignores tf and strata
    cov1 = vv[-c(1L, unlist(sp))] %||% NULL,
    # cov1 = vv[-c(1L, unlist(sp$tf))],
    cov2 = if (is.null(sp$tf))
      # NULL else all.vars(reformulate(vv[sp$tf])),
      NULL else strata(vv[sp$tf], 'tf'),
    
    strata = strata
  )
}

sterms.inner <- function(x) {
  ## survival:::terms.inner
  if (inherits(x, 'formula')) {
    if (length(x) == 3L)
      c(sterms.inner(x[[2L]]), sterms.inner(x[[3L]]))
    else sterms.inner(x[[2L]])
  } else if (class(x) == 'call' &&
             (x[[1L]] != as.name('$') &&
              x[[1L]] != as.name('['))) {
    if (x[[1L]] == '+' || x[[1L]] == '*' || x[[1L]] == '-') {
      if (length(x) == 3L)
        c(sterms.inner(x[[2L]]), sterms.inner(x[[3L]]))
      else sterms.inner(x[[2L]])
    }
    else if (x[[1L]] == as.name('Surv'))
      unlist(lapply(x[-1L], sterms.inner))
    else sterms.inner(x[[2L]])
  } else deparse1(x)
}

strata <- function(formula, pattern = 'strata') {
  ## extract strata from formula
  # cmprsk2:::strata(y ~ x)
  # cmprsk2:::strata(y ~ strata(x) + .)
  # cmprsk2:::strata(y ~ a + b + strata  (factor(x)) + strata( z ) + c)
  # cmprsk2:::strata(y ~ a + b + strata  (factor(x)) + strata(strata( z )) + c)
  if (!inherits(formula, 'formula'))
    formula <- reformulate(formula)
  
  ## remove dot
  dp <- deparse1(formula)
  if (!grepl('strata\\s*\\(', dp))
    return(NULL)
  formula <- as.formula(gsub('\\.\\s+\\+|\\+\\s+\\.', '', dp))
  
  tt <- terms(formula, specials = pattern)
  ii <- attr(tt, 'specials')[[pattern]] - 1L
  tt <- colnames(attr(tt, 'factors'))[ii]
  
  if (length(tt))
    all.vars(reformulate(tt)) else NULL
}

terms.inner <- function(x, survival = FALSE) {
  ## survival:::terms.inner with modifications for crr2/cuminc2 formulae
  if (inherits(x, 'formula')) {
    if (length(x) == 3L) {
      if (!(is.crr2(x) | is.cuminc2(x))) {
        cc <- list(NULL)
        terms2 <- sterms.inner(x[[2L]])
        terms3 <- sterms.inner(x[[3L]])
      } else {
        cc <- strsplit(gsub('^.*\\(|\\)', '', deparse1(x[[2L]])), '==|%in%')
        cc <- trimwsq(cc[[1L]])
        x <- as.formula(sub('\\([^(]*?\\)', '', deparse1(x)))
        terms2 <- terms.inner(x[[2L]])
        terms3 <- sterms.inner(x[[3L]])
      }
      structure(
        list(c(terms2, terms3), cc), dots = any(terms3 == '.')
      )
    } else terms.inner(x[[2L]])
  } else if (class(x) == 'call' &&
             (x[[1L]] != as.name('$') &&
              x[[1L]] != as.name('['))) {
    if (x[[1L]] == '+' || x[[1L]] == '*' || x[[1L]] == '-') {
      if (length(x) == 3L)
        c(terms.inner(x[[2L]]), terms.inner(x[[3L]]))
      else terms.inner(x[[2L]])
    } else if (x[[1L]] == as.name('Surv'))
      unlist(lapply(x[-1L], terms.inner))
    else terms.inner(x[[2L]])
  } else deparse1(x)
}
