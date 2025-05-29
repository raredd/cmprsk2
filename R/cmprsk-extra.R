### cmprsk extra
# lines.cuminc
###


#' Add lines to \code{cuminc}
#' 
#' @param x an object of class \code{\link[cmprsk]{cuminc}}
#' @param lty,col,lwd line type, color, and line width for each curve
#' @param ... additional arguments passed to \code{\link{lines}}
#' 
#' @seealso \code{\link[cmprsk]{plot.cuminc}}
#' 
#' @examples
#' x <- with(transplant, cuminc(futime, event))
#' plot(x)
#' lines(x, col = 'red')
#' 
#' @export

lines.cuminc <- function(x, lty = seq_along(x), col = 1L, lwd = par('lwd'),
                         ...) {
  stopifnot(inherits(x, 'cuminc'))
  nc <- length(x)
  lty <- rep_len(lty, nc)
  col <- rep_len(col, nc)
  lwd <- rep_len(lwd, nc)
  for (ii in seq.int(nc)) {
    lines(x[[ii]][[1L]], x[[ii]][[2L]], lty = lty[ii], col = col[ii], 
          lwd = lwd[ii], ...)
  }
  invisible(x)
}
