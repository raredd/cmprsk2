### utils
# %inside%, %ni%, %||%, islist, assert_class, is.loaded, terms.inner, trimwsq
###


`%inside%` <- function(x, interval) {
  interval <- sort(interval)
  x >= interval[1L] & x <= interval[2L]
}

`%ni%` <- function(x, table) !(match(x, table, nomatch = 0L) > 0L)

`%||%` <- function(x, y) if (is.null(x)) y else x

assert_class <- function(x, class, which = FALSE, message = NULL, warn = FALSE) {
  FUN <- if (warn)
    function(...) warning(..., call. = FALSE)
  else function(...) stop(..., call. = FALSE)
  
  if (is.null(message))
    message <- paste('Object is not of class', toString(shQuote(class)))
  
  if (!all(inherits(x, class, which)))
    FUN(message)
  invisible(TRUE)
}

is.loaded <- function(package) {
  any(grepl(sprintf('package:%s', package), search()))
}

islist <- function(x) inherits(x, 'list')

terms.inner <- function(x, survival = FALSE) {
  ## survival:::terms.inner with modifications
  if (inherits(x, 'formula')) {
    if (length(x) == 3L) {
      cc <- strsplit(gsub('^.*\\(|\\)', '', deparse(x[[2L]])), '==|%in%')[[1L]]
      x <- as.formula(gsub('\\([^(]*?\\)', '', deparse(x)))
      terms3 <- terms.inner(x[[3L]])
      structure(list(c(terms.inner(x[[2L]]), terms3), trimwsq(cc)),
                dots = any(terms3 == '.'))
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
