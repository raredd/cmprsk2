### utils
# %||%, %inside%, assert_class, catlist, islist, is.loaded, nth,
# pvalr, color_pval, rescaler, roundr, signif2, trimwsq
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

color_pval <- function(pv, breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1),
                       cols = colorRampPalette(2:1)(length(breaks)),
                       sig.limit = 0.001, digits = 2L, show.p = FALSE,
                       format_pval = TRUE, na = '-', ...) {
  ## rawr::color_pval -- format/color p-values for printing
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

is.loaded <- function(package) {
  ## accepts strings or ?name (base::is.loaded does not)
  ## otherwise, cannot remember the use of this
  package <- if (is.character(substitute(package)))
    package else deparse1(substitute(package))
  any(grepl(sprintf('package:%s', package), search()))
}

islist <- function(x) {
  ## more explicit is.list (is.list(data.frame()) returns TRUE)
  inherits(x, 'list')
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

trimwsq <- function(x) {
  ## trim white space and quotes
  gsub('^[\'\" \t\r\n]*|[\'\" \t\r\n]*$', '', x)
}
