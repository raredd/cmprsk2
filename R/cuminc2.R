### formula method for cuminc
# cuminc2, plot.cuminc2, print.cuminc2, summary.cuminc2
# 
# unexported:
# get_events, gy_pval, gy_text
###


#' Cumulative incidence analysis
#' 
#' Estimate cumulative incidence functions from competing risks data and test
#' equality across groups.
#' 
#' @param formula a \code{\link[=Surv]{survival object}} formula,
#' \code{Surv(time, status(censor)) ~ response}, where \code{censor} is a
#' unique value of the \code{status} variable indicating the censoring code;
#' note all other unique values of \code{status} will be treated as competing
#' risks
#' @param data a data frame in which to interpret the variables named in
#' \code{formula}
#' @param rho,cencode,subset,na.action passed to \code{\link[cmprsk]{cuminc}};
#' the censoring indicator will be guessed from \code{formula} but may be
#' overridden by \code{cencode}
#' 
#' @examples
#' tp <- within(transplant, {
#'   event_ind <- as.integer(event) - 1L
#'   year <- NULL
#' })
#' 
#' ## no group, no strata
#' form <- Surv(futime, event_ind(0)) ~ 1
#' identical(
#'   cuminc2(form, tp)$cuminc,
#'   with(tp, cuminc(futime, event_ind))
#' )
#' 
#' ## group-only
#' form <- Surv(futime, event_ind(0)) ~ sex
#' identical(
#'   cuminc2(form, tp)$cuminc,
#'   with(tp, cuminc(futime, event_ind, sex))
#' )
#' 
#' ## strata-only
#' form <- Surv(futime, event_ind(0) == 1) ~ strata(abo)
#' identical(
#'   cuminc2(form, tp)$cuminc,
#'   with(tp, cuminc(futime, event_ind, strata = sex))
#' )
#' 
#' ## group and strata
#' form <- Surv(futime, event_ind(0) == 1) ~ sex + strata(abo)
#' identical(
#'   cuminc2(form, tp)$cuminc,
#'   with(tp, cuminc(futime, event_ind, sex, abo))
#' )
#' 
#' 
#' @export

cuminc2 <- function(formula, data, rho = 0, cencode = NULL,
                    subset = NULL, na.action = getOption('na.action')) {
  if (!is.cuminc2(formula))
    stop(
      'Invalid formula - use the following structure or see ?cuminc2:\n',
      '\tSurv(time, status(censor)) ~ response',
      call. = FALSE
    )
  
  form <- parse_formula(formula, data)
  name <- substitute(data)
  Name <- if (length(name) > 1L)
    as.list(name)[[2L]] else name
  na.action <- match.fun(na.action)
  
  subset <- subset %||% rep_len(TRUE, nrow(data))
  data   <- data[subset, ]
  
  idx <- !complete.cases(data[, c(form$lhs, form$rhs)])
  if (any(!!(n <- as.integer(sum(idx))))) {
    message(n, ' observations removed due to missingness', domain = NA)
    data   <- data[!idx, ]
    subset <- subset[!idx]
    assign(as.character(Name), data)
  }
  nr <- nrow(data)
  
  group  <- if (!length(form$rhs))
    rep_len(1L, nr) else do.call('interaction', data[, form$rhs, drop = FALSE])
  strata <- if (is.null(form$strata))
    rep_len(1L, nr) else data[, form$strata]
  
  ci <- cuminc(data[, form$ftime], data[, form$fstatus],
               as.factor(group), as.factor(strata),
               rho, cencode %||% form$cencode, subset, na.action)
  ci <- list(
    cuminc  = ci,
    cuminc2 = setNames(
      cbind.data.frame(
        data[, c(form$ftime, form$fstatus)], group, strata, form$cencode
      ),
      c('time', 'status', 'group', 'strata', 'cencode')
    )
  )
  
  structure(
    ci,
    class = c('cuminc2', 'cuminc')
  )
}

#' \code{cuminc2} plotting method
#' 
#' Plot a \code{\link{cuminc2}} object optional cumulative events or estimate
#' table, Gray's test results, and other features.
#' 
#' @param x an object of class \code{\link{cuminc2}}
#' @param main title of plot
#' @param xlab,ylab x- and y-axis labels
#' @param ... additional parameters (\code{font}, \code{mfrow}, \code{bty},
#' \code{tcl}, \code{cex.lab}, \code{xaxs}, etc) passed to \code{\link{par}}
#' 
#' @seealso
#' \code{\link{cuminc2}}; \code{\link{summary.cuminc2}}
#' 
#' @examples
#' tp <- within(transplant, {
#'   age50 <- factor(+(age > 50))
#'   age_cat <- cut(age, c(0, 40, 60, Inf), c('<40', '40-60', '60+'))
#' })
#' 
#' ci1 <- cuminc2(Surv(futime, event(censored) == ltx) ~ age_cat, tp)
#' plot(ci1)
#' plot(ci1$cuminc)
#' 
#' ci1 <- cuminc2(Surv(futime, event(censored) == death) ~ age50, tp)
#' plot(ci1, lty.ci = c(1,1,2,2,3,3), col.ci = 1:2)
#' 
#' ci2 <- cuminc2(Surv(futime, event(censored) == death) ~ 1, tp)
#' plot(ci2, wh.events = 'est')
#' 
#' @export

plot.cuminc2 <- function(x,
                         main = NULL,
                         xlab = 'Time', ylab = 'Probability',
                         labels = names(xx),
                         strata.lab = labels,
                         xlim = NULL, ylim = NULL,
                         col.ci = seq_along(xx),
                         lty.ci = par('lty'), lwd.ci = par('lwd'),
                         cex.axis = par('cex.axis'),
                         gy_test = TRUE,
                         events = TRUE, col.events = col.ci,
                         wh.events = c('events', 'est', 'est.sd', 'est.ci'),
                         events.digits = 3L,
                         events.lab = NULL,
                         events.lines = TRUE, events.col = FALSE,
                         xaxis.at = pretty(xlim),
                         yaxis.at = pretty(ylim),
                         xaxis.lab = xaxis.at, yaxis.lab = yaxis.at,
                         events.at = xaxis.at,
                         strata.order = seq_along(xx),
                         extra.margin = 5L,
                         panel.first = NULL, panel.last = NULL,
                         ...) {
  assert_class(x, 'cuminc2')
  
  ci <- x[['cuminc']]
  xx <- ci[names(ci) != 'Tests']
  ng <- length(xx)
  ng <- max(ng, 1L)
  
  col.ci <- rep_len(col.ci, ng)
  lwd.ci <- rep_len(lwd.ci, ng)
  lty.ci <- rep_len(lty.ci, ng)
  
  wh.events <- match.arg(wh.events)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  par(mar = c(4 + ng * events,
              4 + pmax(4, extra.margin) - 3 * !events,
              2,
              2 + 3 * events))
  
  time <- unlist(sapply(xx, '[[', 'time'))
  if (is.null(xlim))
    xlim <- c(0, max(time))
  if (is.null(ylim))
    ylim <- c(0, 1)
  
  plot(xx[[1L]][[1L]], xx[[1L]][[2L]], type = 'n',
       xlim = xlim, ylim = ylim, ann = FALSE, axes = FALSE,
       panel.first = panel.first,
       panel.last = {
         box(bty = par('bty'))
         axis(1L, xaxis.at, FALSE, lwd = 0, lwd.ticks = 1)
         axis(1L, xaxis.at, xaxis.lab, FALSE, -0.5, cex.axis = cex.axis)
         axis(2L, yaxis.at, yaxis.lab, las = 1L, cex.axis = cex.axis)
         title(xlab = xlab, line = 1.5, adj = 0.5)
         title(ylab = ylab, main = main)
       })
  
  ## ci lines
  for (ii in seq.int(ng))
    lines(xx[[ii]][[1L]], xx[[ii]][[2L]], ...,
          lty = lty.ci[ii], col = col.ci[ii], lwd = lwd.ci[ii])
  
  if (events) {
    usr <- par('usr')
    
    events.at <- events.at[events.at < usr[2L]]
    
    ## set colors for lines of text
    col.events <- if (isTRUE(events.col))
      col.ci else if (!identical(events.col, FALSE) &&
                      length(events.col) == ng)
        events.col else rep_len(1L, ng)
    
    ## at-risk label
    if (!(identical(events.lab, FALSE))) {
      events.lab <- if (is.null(events.lab))
        switch(wh.events,
               events = 'Cumulative events',
               est    = 'Estimate',
               est.sd = 'Estimate +/- Std. dev',
               est.ci = 'Estimate [LCI, UCI]')
      else events.lab
      
      mtext(events.lab, side = 1L, at = usr[1L], line = 1.5,
            adj = 1, col = 1L, las = 1L, cex = cex.axis, font = 2L)
    }
    
    ## labels for each row in at-risk table
    group.name.pos <- diff(usr[1:2]) / -8
    padding  <- abs(group.name.pos / 8)
    line.pos <- seq.int(ng)[order(strata.order)] + 2L
    
    mtext(strata.lab, side = 1L, line = line.pos, adj = 1, las = 1L,
          col = col.events, at = group.name.pos, cex = cex.axis, font = 2L)
    
    ## labels for total events
    total.lab <- 'Total events'
    n.events <- x$cuminc2
    n.events <- n.events[n.events$status %ni% n.events$cencode, ]
    n.events <- with(droplevels(n.events), table(group, status))
    mtext(c(n.events), side = 1L, at = abs(group.name.pos) + usr[2L],
          line = line.pos, adj = 1, col = col.events, las = 1L,
          cex = cex.axis, font = 2L)
    mtext(total.lab, side = 1L, at = abs(group.name.pos) + usr[2L],
          line = 1.5, adj = 1, col = 1L, las = 1L, cex = cex.axis,
          font = 2L)
    
    ## draw matching lines for n events  
    if (events.lines)
      for (ii in seq.int(ng))
        ## mess with the 4 here to adjust the length of the events.line
        axis(1L, c(group.name.pos + padding, 0 - 4 * padding), xpd = NA,
             labels = FALSE, line = line.pos[ii] + 0.6, lwd.ticks = 0,
             col = col.ci[ii], lty = lty.ci[ii], lwd = lwd.ci[ii])
    
    ## number of events
    ss <- switch(
      wh.events,
      events = summary(x, times = events.at)$events,
      est    = timepoints2(x, times = events.at, sd = FALSE, html = FALSE,
                           ci = FALSE, digits = events.digits),
      est.sd = timepoints2(x, times = events.at, sd = TRUE, html = FALSE,
                           ci = FALSE, digits = events.digits),
      est.ci = timepoints2(x, times = events.at, sd = FALSE, html = FALSE,
                           ci = TRUE, digits = events.digits)
    )
    d1 <- data.frame(time = rep(as.numeric(colnames(ss)), each = nrow(ss)),
                     n.risk = c(ss), strata = rep(rownames(ss), ncol(ss)),
                     stringsAsFactors = FALSE)
    d1$strata <- factor(d1$strata, rownames(ss))
    # d1 <- d1[d1$n.risk > 0, ]
    d2 <- split(d1, d1$strata)
    
    ## right-justify numbers
    ndigits <- lapply(d2, function(x) nchar(x[, 2L]))
    max.len <- max(sapply(ndigits, length))
    L <- do.call('rbind', lapply(ndigits, `length<-`, max.len))
    nd <- apply(L, 2L, max, na.rm = TRUE)
    
    for (ii in seq.int(ng)) {
      tmp <- d2[[ii]]
      w.adj <- strwidth('0', cex = cex.axis, font = par('font')) /
        2 * nd[seq.int(nrow(tmp))]
      mtext(tmp$n.risk, side = 1L, at = tmp$time + w.adj, cex = cex.axis,
            las = 1L, line = line.pos[ii], adj = 1, col = col.events[ii])
    }
  }
  
  ## add gray text in upper right corner
  if (gy_test) {
    txt <- tryCatch(gy_text(ci), error = function(e) 'n/a')
    if (is.null(txt))
      invisible(NULL)
    else legend('topleft', legend = paste0(names(txt), ': ', txt),
                bty = 'n')
  }
  
  panel.last
  
  invisible(NULL)
}

#' \code{cuminc2} print method
#' 
#' @param x an object of class \code{\link{cuminc2}}
#' @param ... ignored
#' 
#' @seealso
#' \code{\link{cuminc2}}; \code{\link{summary.cuminc2}}
#' 
#' @export

print.cuminc2 <- function(x, ...) {
  if (inherits(x, 'cuminc2')) {
    print(x[['cuminc']])
  } else print(x)
  invisible(x)
}

#' \code{cuminc2} summary method
#'
#' @param object an object of class \code{\link{cuminc2}}
#' @param times a vector of times
#' @param ... ignored
#' 
#' @seealso
#' \code{\link{cuminc2}}; \code{\link{print.cuminc2}};
#' \code{\link[cmprsk]{timepoints}}; \code{\link{timepoints2}}
#' 
#' @return
#' A list with the following components:
#' 
#' \item{\code{est}}{a matrix of estimates of the subdistributions for the
#' cumulative incidence model(s); see \code{\link[cmprsk]{timepoints}}}
#' \item{\code{var}}{a matrix giving the corresponding variances for the
#' cumulative incidence model(s); see \code{\link[cmprsk]{timepoints}}}
#' \item{\code{events}}{a matrix with the number of cumulative events up
#' to and including each value in \code{times}}
#' \item{\code{total_events}}{a matrix giving the total number of events}
#' 
#' @examples
#' ci <- cuminc2(Surv(futime, event(censored)) ~ sex + strata(abo),
#'               data = transplant)
#' summary(ci)
#' summary(ci, times = c(0, 250, 500))
#'
#' @export

summary.cuminc2 <- function(object, times = NULL, ...) {
  assert_class(object, 'cuminc2')
  
  x <- object$cuminc2
  x <- droplevels(x[x$status %ni% x$cencode, ])
  
  st <- levels(as.factor(x$group))
  
  sp <- split(x, interaction(x$status, x$strata, drop = TRUE))
  sp <- split(x, interaction(x$status, drop = TRUE))
  
  sp <- lapply(sp, function(x) {
    x <- x[order(x$time), ]
    x[, st] <- lapply(st, function(y) {
      cumsum(x$group %in% y)
    })
    x
  })
  
  if (is.null(times))
    times <- pretty(x$time)
  tp <- object[['cuminc']]
  tp <- timepoints(tp, times)
  
  times <- colnames(tp$est)
  sp <- lapply(sp, function(x) {
    get_events(x[, st, drop = FALSE], x$time, times)
  })
  
  res <- do.call('rbind', sp)
  rownames(res) <- paste(st, rep(names(sp), each = length(st)))
  colnames(res) <- times
  total_events <- matrix(table(x$group, x$status),
                         dimnames = list(rownames(res), 'Events'))
  ## "atrisk"
  # res <- c(total_events) - res
  
  c(tp, list(events = res, total_events = total_events))
}

get_events <- function(data, time, timepoints) {
  timepoints <- as.numeric(timepoints)
  
  res <- lapply(timepoints, function(x) {
    idx <- time <= x
    res <- if (x < 0)
      matrix(0L, 1L, ncol(data))
    else if (sum(idx))
      data[max.col(t(idx), ties.method = 'last'), ]
    else data[1L, ]
    as.matrix(res)
  })
  
  t(do.call('rbind', res))
}

gy_pval <- function(x, details = FALSE) {
  # gy_pval(cuminc2); gy_pval(cuminc2$cuminc)
  assert_class(x, c('cuminc2', 'cuminc'))
  x <- x[['cuminc']] %||% x
  stopifnot(
    inherits(x, 'cuminc')
  )
  if (details)
    x$Tests else x$Tests[, 2L]
}

gy_text <- function(x, ...) {
  # gy_text(cuminc2)
  if (inherits(x, c('cuminc2', 'cuminc')))
    x <- gy_pval(x, TRUE)
  
  if (is.null(x))
    return(invisible(NULL))
  
  txt <- apply(x, 1L, function(x) {
    sprintf('%s (%s df), %s', sprintf('%.2f', x['stat']),
            x['df'], pvalr(x['pv'], show.p = TRUE))
  })
  
  # bquote(paste(chi^2, ' = ', .(txt)))
  txt
}
