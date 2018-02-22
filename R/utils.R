### utils
# %||%, %inside%, %ni%, assert_class, islist, is.loaded, nth, parse_formula,
# pvalr, strata, terms.inner, sterms.inner, trimwsq
# 
# formula s3 methods:
# is.crr2, is.crr2.default, is.crr2.formula, is.cuminc2, is.cuminc2.default,
# is.cuminc2.formula
###


`%||%` <- function(x, y) {
  if (is.null(x) | !length(x)) y else x
}

`%inside%` <- function(x, interval) {
  interval <- sort(interval)
  x >= interval[1L] & x <= interval[2L]
}

`%ni%` <- function(x, table) {
  !(match(x, table, nomatch = 0L) > 0L)
}

assert_class <- function(x, class, which = FALSE,
                         message = NULL, warn = FALSE) {
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
  paste0(paste(names(x), x, sep = ' = ', collapse = ', '))
}

is.crr2 <- function(x) {
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
  x <- paste(x, collapse = ' ')
  
  ## checks for status(0) but not status(0) == 1
  # grepl('[^(]+\\([^(~]+\\((?=.+~|[^~]+$)', x, perl = TRUE)
  
  ## asserts status(0) == 1
  grepl('Surv\\([^(]+\\([^(]+\\)\\s*(==|%in%)[^)]+\\)', x)
}

is.cuminc2 <- function(x) {
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
  
  ## status(0) == 1 equality is optional
  grepl('Surv\\([^(]+\\([^(]+\\)\\s*(?:(==|%in%)[^)]+)?\\)', x)
}

islist <- function(x) {
  ## is.list(data.frame()) returns TRUE
  inherits(x, 'list')
}

is.loaded <- function(package) {
  package <- if (is.character(substitute(package)))
    package else deparse(substitute(package))
  any(grepl(sprintf('package:%s', package), search()))
}

parse_formula <- function(formula, data = NULL) {
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
  
  ftime   <- lhs[1L]
  fstatus <- lhs[2L]
  
  cencode  <- form[[2L]][1L] %||% NULL
  failcode <- if (is.na(failcode <- form[[2L]][2L]))
    NULL else failcode
  
  if (!is.null(data)) {
    dname <- deparse(substitute(data))
    vars <- c(lhs, rhs, strata)
    if (length(vars <- vars[vars %ni% names(data)]))
      stop(
        sprintf('%s not found in %s',
                toString(shQuote(vars)), shQuote(dname))
      )
    
    if (!is.numeric(na.omit(data[, ftime])) |
        any(na.omit(data[, ftime]) < 0))
      stop(
        sprintf('%s should be numeric values > 0', shQuote(ftime))
      )
    
    codes <- c(cencode, failcode)
    if (length(codes <- codes[codes %ni% data[, fstatus]]))
      stop(
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

## rawr::pvalr
pvalr <- function(pvals, sig.limit = 0.001, digits = 3L,
         html = FALSE, show.p = FALSE) {
  stopifnot(sig.limit > 0, sig.limit < 1)
  show.p <- show.p + 1L
  html   <- html + 1L
  
  sapply(pvals, function(x, sig.limit) {
    if (is.na(x) | !nzchar(x))
      return(NA)
    if (x >= 0.99)
      return(paste0(c('','p ')[show.p], c('> ','&gt; ')[html], '0.99'))
    if (x >= 0.9)
      return(paste0(c('','p ')[show.p], c('> ','&gt; ')[html], '0.9'))
    if (x < sig.limit) {
      paste0(c('', 'p ')[show.p], c('< ', '&lt; ')[html], format(sig.limit))
    } else {
      nd <- c(digits, 2L, 1L)[findInterval(x, c(-Inf, .1, .5, Inf))]
      paste0(c('','p = ')[show.p], roundr(x, nd))
    }
  }, sig.limit)
}

nth <- function(x, p, n = NULL, keep_split = FALSE, repl = '$$$', ...) {
  # s <- 'this  is a  test string to use   for testing   purposes'
  # nth(s, '\\s+')
  # nth(s, '\\s+', 3)
  # nth(s, '\\s+', c(3, 5))
  # nth(s, '\\s+', c(3, 5), keep_split = TRUE)
  # nth(s, '\\s{2,}', keep_split = TRUE)
  # nth(s, '\\s{2,}', 2:4)
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
  (x - from[1L]) / diff(from) * diff(to) + to[1L]
}

roundr <- function(x, digits = 1L) {
  sprintf('%.*f', digits, x)
}

strata <- function(formula) {
  formula <- deparse(formula)
  strata  <- trimws(gsub('strata\\s*\\(([^)]+)|.', '\\1', formula))
  if (nzchar(strata))
    strata else NULL
}

## rawr::tcol
tcol <- function(colors, trans = NULL, alpha = NULL) {
  trans <- trans %||% 255L
  stopifnot(
    trans %inside% c(0L, 255L) | is.na(trans)
  )
  
  ## convert alpha to trans
  if (!is.null(alpha)) {
    stopifnot(
      alpha %inside% c(0, 1) | is.na(alpha)
    )
    trans <- as.integer(rescaler(alpha, to = c(0, 255), from = c(0, 1)))
  }
  
  ## get color and trans to conformable lengths
  if (length(colors) != length(trans) & 
      !any(c(length(colors), length(trans)) == 1L))
    stop('Vector lengths are not conformable')
  if (length(colors) == 1L & length(trans) > 1L)
    colors <- rep_len(colors, length(trans))
  if (length(trans) == 1L & length(colors) > 1L)
    trans <- rep_len(trans, length(colors))
  
  ## if color == 0, returns NA
  if (length(nocol <- which(colors == 0))) {
    colors[nocol] <- 1
    trans[nocol] <- NA
  }
  
  res <- paste0('#', apply(apply(rbind(col2rgb(colors)), 2L, function(x)
    format(as.hexmode(x), width = 2L)), 2L, paste, collapse = ''))
  res <- Map(paste0, res, tryCatch(
    as.character(as.hexmode(trans)),
    error = function(e) '', warning = function(w) ''
  ))
  res <- unname(unlist(res))
  
  ## return NAs and/or set color to transparent
  res[is.na(colors) | is.na(trans)] <- NA
  res[colors %in% 'transparent'] <- 'transparent'
  
  res
}

terms.inner <- function(x, survival = FALSE) {
  ## survival:::terms.inner with modifications
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
  gsub('^[\'\" \t\r\n]*|[\'\" \t\r\n]*$', '', x)
}
