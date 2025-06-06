% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmprsk.R
\name{lines.cuminc}
\alias{lines.cuminc}
\title{Add lines to \code{cuminc}}
\usage{
\method{lines}{cuminc}(x, lty = seq_along(x), col = 1L, lwd = par("lwd"), ...)
}
\arguments{
\item{x}{an object of class \code{\link[cmprsk]{cuminc}}}

\item{lty, col, lwd}{line type, color, and line width for each curve}

\item{...}{additional arguments passed to \code{\link{lines}}}
}
\description{
Add lines to \code{cuminc}
}
\examples{
x <- with(transplant, cuminc(futime, event))
plot(x)
lines(x, col = 'red')

}
\seealso{
\code{\link[cmprsk]{plot.cuminc}}
}
