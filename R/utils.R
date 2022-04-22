### utils
# %||%, %inside%, %ni%, assert_class, catlist, islist, is.loaded, nth,
# parse_formula, pvalr, color_pval, rescaler, roundr, signif2, strata,
# terms.inner, sterms.inner, trimwsq
# 
# formula s3 methods:
# is.crr2, is.crr2.default, is.crr2.formula, is.cuminc2, is.cuminc2.default,
# is.cuminc2.formula
###


`%||%` <- function(x, y) {
  ## return y if x is null
  if (is.null(x) || !length(x))
    y else x
}

`%inside%` <- function(x, interval) {
  ## check if x is weakly inside interval
  interval <- sort(interval)
  x >= interval[1L] & x <= interval[2L]
}

`%ni%` <- function(x, table) {
  ## negation of %in%
  !(match(x, table, nomatch = 0L) > 0L)
}

assert_class <- function(x, class, which = FALSE,
                         message = NULL, warn = FALSE) {
  ## checking classes
  name <- substitute(x)
  FUN <- if (warn)
    warning else stop
  formals(FUN)$call. <- FALSE
  
  if (is.null(message))
    message <- paste(shQuote(name), 'is not of class',
                     toString(shQuote(class)))
  
  if (!all(inherits(x, class, which)))
    FUN(message)
  
  invisible(TRUE)
}

catlist <- function(x) {
  ## concatenate a list to print
  paste0(paste(names(x), x, sep = ' = ', collapse = ', '))
}

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
    deparse(x) else deparse(x[[2L]])
  x <- paste(x, collapse = '')
  
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
    deparse(x) else deparse(x[[2L]])
  x <- paste(x, collapse = '')
  
  ## status(0) == 1 equality is optional
  grepl('Surv\\([^(]+\\([^(]+\\)\\s*(?:(==|%in%)[^)]+)?\\)', x) &
    grepl('~', x)
}

islist <- function(x) {
  ## more explicit is.list (is.list(data.frame()) returns TRUE)
  inherits(x, 'list')
}

is.loaded <- function(package) {
  ## accepts strings or ?name (base::is.loaded does not)
  ## otherwise, cannot remember the use of this
  package <- if (is.character(substitute(package)))
    package else deparse(substitute(package))
  any(grepl(sprintf('package:%s', package), search()))
}

parse_formula <- function(formula, data = NULL, name = NULL) {
  ## parse crr2/cuminc2 formulae into useful pieces
  if (!(is.crr2(formula) | is.cuminc2(formula)))
    stop(
      'Invalid formula - see ?crr2 or ?cuminc2 for deatails'
    )
  
  term <- terms.inner(formula)
  form <- rapply(term, trimwsq, how = 'list')
  
  lhs  <- form[[1L]][1:2]
  rhs  <- if (attr(term, 'dots'))
    setdiff(names(data), lhs) else form[[1L]][-(1:2)]
  if (identical(rhs, '1'))
    rhs <- NULL
  
  strata <- strata(formula)
  rhs    <- setdiff(rhs, strata) %||% NULL
  rhs    <- if (is.null(strata))
    rhs else setdiff(rhs, 'strata')
  
  ftime   <- lhs[1L]
  fstatus <- lhs[2L]
  
  cencode  <- form[[2L]][1L] %||% NULL
  failcode <- if (is.na(failcode <- form[[2L]][2L]))
    NULL else failcode
  
  if (!is.null(data)) {
    dname <- name %||% deparse(substitute(data))
    vars <- c(lhs, rhs, strata)
    
    if (length(vars <- vars[vars %ni% names(data)]))
      stop(
        sprintf('%s not found in %s',
                toString(shQuote(vars)), shQuote(dname))
      )
    
    if (!is.numeric(na.omit(data[, ftime, drop = TRUE])) |
        any(na.omit(data[, ftime, drop = TRUE]) < 0))
      stop(
        sprintf('%s should be numeric values >= 0', shQuote(ftime))
      )
    
    codes <- c(cencode, failcode)
    if (length(codes <- codes[codes %ni% data[, fstatus, drop = TRUE]]))
      warning(
        sprintf('%s not found in %s[, %s]',
                toString(shQuote(codes)), dname, shQuote(fstatus))
      )
  }
  
  list(
    lhs = lhs,
    rhs = rhs,
    
    ftime   = ftime,
    fstatus = fstatus,
    
    cencode  = cencode,
    failcode = failcode,
    
    strata = strata
  )
}

pvalr <- function(pv, sig.limit = 0.001, digits = 2L, scientific = FALSE,
                  html = FALSE, show.p = FALSE, ...) {
  ## rawr::pvalr -- format p-values for printing
  stopifnot(
    sig.limit > 0,
    sig.limit < 1
  )
  
  pstr <- c('', 'p ')[show.p + 1L]
  high <- 1 - 1 / 10 ^ digits
  
  sapply(pv, function(x) {
    if (is.na(x) | !nzchar(x) | !is.numeric(x))
      NA
    else if (x > high)
      paste0(pstr, c('> ', '&gt; ')[html + 1L], high)
    else if (x < sig.limit)
      paste0(pstr, c('< ', '&lt; ')[html + 1L],
             format.pval(sig.limit, scientific = scientific))
    else paste0(c('', 'p = ')[show.p + 1L], signif2(x, digits))
  })
}

color_pval <- function(pv, breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1),
                       cols = colorRampPalette(2:1)(length(breaks)),
                       sig.limit = 0.001, digits = 2L, show.p = FALSE,
                       format_pval = TRUE, na = '-', ...) {
  if (!is.numeric(pv))
    return(pv)
  pvn <- pv
  
  stopifnot(length(breaks) == length(cols))
  
  pv <- if (isTRUE(format_pval))
    pvalr(pvn, sig.limit, digits, scientific = FALSE, html = TRUE, show.p)
  else if (identical(format_pval, FALSE))
    pv else format_pval(pv)
  
  pvc <- cols[findInterval(pvn, breaks)]
  res <- sprintf('<font color=\"%s\">%s</font>', pvc, pv)
  
  replace(res, grepl('>NA<', res, fixed = TRUE), na)
}

nth <- function(x, p, n = NULL, keep_split = FALSE, repl = '$$$', ...) {
  ## split a string, s, at the n-th occurrence of a pattern, p
  # s <- 'this  is a  test string to use   for testing   purposes'
  # cmprsk2:::nth(s, '\\s+')
  # cmprsk2:::nth(s, '\\s+', 3)
  # cmprsk2:::nth(s, '\\s{2,}', 3) ## compare
  # cmprsk2:::nth(s, '\\s+', c(3, 5))
  # cmprsk2:::nth(s, '\\s+', c(3, 5), keep_split = TRUE)
  # cmprsk2:::nth(s, '\\s{2,}', keep_split = TRUE)
  # cmprsk2:::nth(s, '\\s{2,}', 2:4)
  stopifnot(
    is.character(x),
    is.character(p),
    length(x) == 1L
  )
  m <- gregexpr(p, x, ...)
  l <- length(attr(m[[1L]], 'match.length'))
  n <- if (is.null(n))
    seq.int(l) else n[n < l]
  
  regmatches(x, m)[[1L]][n] <- if (keep_split)
    paste0(repl, regmatches(x, m)[[1L]][n], repl) else repl
  
  strsplit(x, repl, fixed = TRUE)[[1L]]
}

rescaler <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
  ## stripped version of scales::rescale
  (x - from[1L]) / diff(from) * diff(to) + to[1L]
}

roundr <- function(x, digits = 1L) {
  ## round without dropping trailing 0s (returns character strings)
  sprintf('%.*f', digits, x)
}

signif2 <- function(x, digits = 6L) {
  sapply(x, function(xx) {
    s <- signif(xx, digits = digits)
    formatC(s, digits = digits, format = 'fg', flag = '#')
  })
}

strata <- function(formula) {
  ## extract strata from formula
  # cmprsk2:::strata(y ~ x)
  # cmprsk2:::strata(y ~ strata(x))
  # cmprsk2:::strata(y ~ a + b + strata  (x) + strata( z ) + c)
  formula <- deparse(formula)
  formula <- paste0(formula, collapse = '')
  
  m <- gregexpr('strata\\s*\\(\\K([^)]+)', formula, perl = TRUE)
  strata <- trimws(regmatches(formula, m)[[1L]])
  
  if (length(strata))
    strata else NULL
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
        cc <- strsplit(gsub('^.*\\(|\\)', '', deparse(x[[2L]])), '==|%in%')
        cc <- trimwsq(cc[[1L]])
        x <- as.formula(sub('\\([^(]*?\\)', '', deparse(x)))
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
  } else deparse(x)
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
  } else (deparse(x))
}

trimwsq <- function(x) {
  ## trim white space and quotes
  gsub('^[\'\" \t\r\n]*|[\'\" \t\r\n]*$', '', x)
}
