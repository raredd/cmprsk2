### utils
# %ni%, %||%, islist, assert_class, is.loaded, terms.inner, trimwsq
###


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

terms.inner <- function(x) {
  ## survival:::terms.inner with modifications
  if (inherits(x, 'formula')) {
    if (length(x) == 3L) {
      cc <- strsplit(gsub('^.*\\(|\\)', '', deparse(x[[2]])), '==|%in%')[[1]]
      x <- as.formula(gsub('\\([^(]*?\\)', '', deparse(x)))
      terms3 <- terms.inner(x[[3]])
      structure(list(c(terms.inner(x[[2]]), terms3), trimwsq(cc)),
                dots = any(terms3 == '.'))
    } else terms.inner(x[[2]])
  } else if (class(x) == 'call' &&
             (x[[1]] != as.name('$') &&
              x[[1]] != as.name('['))) {
    if (x[[1]] == '+' || x[[1]] == '*' || x[[1]] == '-') {
      if (length(x) == 3L)
        c(terms.inner(x[[2]]), terms.inner(x[[3]]))
      else terms.inner(x[[2]])
    } else if (x[[1]] == as.name('Surv'))
      unlist(lapply(x[-1], terms.inner))
    else terms.inner(x[[2]])
  } else deparse(x)
}

trimwsq <- function(x) gsub('^[\'\" \t\r\n]*|[\'\" \t\r\n]*$', '', x)
