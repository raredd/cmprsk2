### cmprsk::cuminc functions
# timepoints2, print.timepoints2
###


#' \code{timepoints} formatter
#' 
#' Format output from \code{\link[cmprsk]{timepoints}} into a single matrix
#' with estimate +/- standard deviation.
#' 
#' @param w,times arguments passed to \code{\link{timepoints}}
#' @param digits number of digits past the decimal point to keep
#' @param sd logical; if \code{FALSE}, the standard deviation will not be
#' shown with the estimate
#' @param ci logical; not implemented
#' @param html logical; if \code{TRUE}, an html-friendly format is returned;
#' the print method for \code{timepoints2} will use \code{\link{htmlTable}}
#' if \code{html = TRUE}
#' @param x an object of class \code{"timepoints2"}
#' @param ... additional arguments passed to or from other methods
#' 
#' @examples
#' ci <- cuminc2(Surv(futime, event(censored)) ~ sex, transplant)
#' timepoints2(ci, sd = TRUE)
#' timepoints2(ci, html = TRUE)
#' 
#' print(
#'   timepoints2(ci, sd = TRUE, html = TRUE),
#'   rgroup = c('Death', 'Ltx', 'Withdraw'), n.rgroup = c(2, 2, 2),
#'   caption = 'Timepoints<sup>&dagger;</sup>',
#'   tfoot = '<sup>&dagger;</sup>Estimate &pm; SD'
#' )
#' 
#' 
#' ## example from cmprsk::cuminc
#' set.seed(2)
#' ss <- rexp(100)
#' gg <- factor(sample(1:3,100,replace=TRUE),1:3,c('a','b','c'))
#' cc <- sample(0:2,100,replace=TRUE)
#' strt <- sample(1:2,100,replace=TRUE)
#' xx <- cuminc(ss,cc,gg,strt)
#' 
#' timepoints(xx, times = 0:4)
#' timepoints2(xx, times = 0:4, digits = 5)
#' 
#' @export

timepoints2 <- function(w, times, digits = 3L,
                        sd = FALSE, ci = FALSE, html = FALSE) {
  w <- if (inherits(w, 'cuminc2'))
    w[['cuminc']]
  else if (inherits(w, 'cuminc'))
    w else stop('\'w\' should be of class \'cuminc\' or \'cuminc2\'')
  tt <- na.omit(unlist(sapply(w, `[`, 'time')))
  rr <- range(tt)
  
  if (missing(times))
    times <- pretty(tt)
  times <- times[times %inside% rr]
  tp <- timepoints(w, times)
  res <- tp$est
  
  fmt <- if (ci)
    sprintf('%%.%sf [%%.%sf - %%.%sf]',
            digits, digits, digits)
  else {
    sprintf('%%.%sf %s %%.%sf',
            digits, if (html) '&pm;' else '+/-', digits)
  }
  
  if (ci) {
    res[] <- mapply(function(x, y, z)
      sprintf(fmt, x, y, z),
      tp$est,
      pmax(tp$est - 1.96 * sqrt(tp$var), 0),
      tp$est + 1.96 * sqrt(tp$var))
  } else {
    res[] <- mapply(function(x, y)
      sprintf(fmt, x, y), tp$est, sqrt(tp$var))
    if (!sd)
      res <- gsub(' [+&].*$', '', res)
  }
  
  res <- gsub('NA.*$', '-', res)
  
  attr(res, 'html') <- html
  
  structure(
    res,
    class = 'timepoints2'
  )
}

#' @rdname timepoints2
#' @export
print.timepoints2 <- function(x, ..., html = attr(x, 'html')) {
  assert_class(x, 'timepoints2')
  x <- if (html) {
    structure(
      htmlTable::htmlTable(
        x, css.cell = 'padding: 0px 5px 0px; white-space: nowrap;', ...),
      class = 'htmlTable'
    )
  } else {
    class(x) <- NULL
    attr(x, 'html') <- NULL
    x
  }
  
  print(x)
}
