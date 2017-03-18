### utils
# %inside%, %ni%, %||%, islist, assert_class, is.loaded, terms.inner, trimwsq
###


`%inside%` <- function(x, interval) {
  interval <- sort(interval)
  x >= interval[1L] & x <= interval[2L]
}

`%ni%` <- function(x, table) !(match(x, table, nomatch = 0L) > 0L)

`%||%` <- function(x, y) if (is.null(x)) y else x

assert_class <- function(x, class, which = FALSE,
                         message = NULL, warn = FALSE) {
  name <- substitute(x)
  FUN <- if (warn)
    function(...) warning(..., call. = FALSE)
  else function(...) stop(..., call. = FALSE)
  
  if (is.null(message))
    message <- paste(shQuote(name), 'is not of class',
                     toString(shQuote(class)))
  
  if (!all(inherits(x, class, which)))
    FUN(message)
  invisible(TRUE)
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
  
  ## checks for status(0) but not status(0) == 1
  # grepl('[^(]+\\([^(~]+\\((?=.+~|[^~]+$)', x, perl = TRUE)
  
  ## asserts status(0) == 1
  grepl('Surv\\([^(]+\\([^(]+\\)\\s*(==|%in%)[^)]+\\)', x)
}

is.loaded <- function(package) {
  any(grepl(sprintf('package:%s', package), search()))
}

islist <- function(x) inherits(x, 'list')

terms.inner <- function(x, survival = FALSE) {
  ## survival:::terms.inner with modifications
  if (inherits(x, 'formula')) {
    if (length(x) == 3L) {
      if (!is.crr2(x)) {
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
      structure(list(c(terms2, terms3), cc), dots = any(terms3 == '.'))
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

trimwsq <- function(x) gsub('^[\'\" \t\r\n]*|[\'\" \t\r\n]*$', '', x)
