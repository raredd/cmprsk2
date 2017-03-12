### cuminc functions
# timepoints2, print.timepoints2
###


#' \code{timepoints} formatter
#' 
#' Format output from \code{\link[cmprsk]{timepoints}} into a single matrix
#' with estimate +/- standard deviation.
#' 
#' @param w,times arguments passed to \code{\link{timepoints}}
#' @param digits number of digits past the decimal point to keep
#' @param ci logical; not implemented
#' @param html logical; if \code{TRUE}, html-friendly formatted is returned;
#' the print method for \code{timepoints2} will use \code{\link{htmlTable}}
#' if \code{html = TRUE}
#' @param x an object of class \code{"timepoints2"}
#' @param ... additional arguments passed to or from other methods
#' 
#' @examples
#' ## example from cmprsk::cuminc
#' set.seed(2)
#' ss <- rexp(100)
#' gg <- factor(sample(1:3,100,replace=TRUE),1:3,c('a','b','c'))
#' cc <- sample(0:2,100,replace=TRUE)
#' strt <- sample(1:2,100,replace=TRUE)
#' print(xx <- cuminc(ss,cc,gg,strt))
#' 
#' timepoints(xx, times = 0:4)
#' timepoints2(xx, html = FALSE)
#' 
#' print(timepoints2(xx),
#'       rgroup = c('One', 'Two'), n.rgroup = c(3,3),
#'       caption = 'Timepoints<sup>&dagger;</sup>',
#'       tfoot = '<sup>&dagger;</sup>Est &pm; sd')
#' 
#' @export

timepoints2 <- function(w, times, digits = 3L, ci = TRUE, html = TRUE) {
  stopifnot(inherits(w, 'cuminc'))
  tt <- na.omit(unlist(sapply(w, `[`, 'time')))
  rr <- range(tt)
  
  if (missing(times))
    times <- pretty(tt)
  times <- times[times %inside% rr]
  tp <- timepoints(w, times)
  res <- tp$est
  
  fmt <- sprintf('%%.~f %s %%.~f', if (html) '&pm;' else '+/-')
  fmt <- gsub('~', digits, fmt)
  res[] <- mapply(function(x, y) sprintf(fmt, x, y), tp$est, sqrt(tp$var))
  
  attr(res, 'html') <- html
  structure(res, class = 'timepoints2')
}

#' @rdname timepoints2
#' @export
print.timepoints2 <- function(x, ..., html = attr(x, 'html')) {
  stopifnot(inherits(x, 'timepoints2'))
  x <- if (html) {
    x <- htmlTable::htmlTable(
      x, css.cell = 'padding: 0px 5px 0px; white-space: nowrap;', ...
    )
    ## bug in htmlTable v1.9 with class == c('html', ...)
    class(x) <- 'htmlTable'
    x
  } else {
    class(x) <- NULL
    attr(x, 'html') <- NULL
    x
  }
  print(x)
}
