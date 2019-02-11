### plotting functions
# ciplot, plot.cuminc2, ciplot_by
###


#' \code{cuminc2} plotting method
#' 
#' Plot a \code{\link{cuminc2}} object with optional cumulative events or
#' estimate table, Gray's test results, and other features.
#' 
#' @param x an object of class \code{\link{cuminc2}}
#' @param col.ci,lty.ci,lwd.ci line color, type, and width for each curve
#' @param events logical; if \code{TRUE}, a cumulative events table is drawn
#' @param atrisk logical; if \code{TRUE}, the number at risk is added to the
#' events table
#' @param events.total logical or numeric; if \code{TRUE}, cumulative events
#' for each group is added to events table at a calculated position; for more
#' control, use a specific x-coordinate
#' @param wh.events a character string giving the type of \code{events} table
#' to show; one of \code{"events"} (cumulative number of events), \code{"est"}
#' (estimates, see \code{\link[cmprsk]{timepoints}}), \code{"est.sd"} (estimate
#' +/- standard deviation), or \code{"est.ci"} (estimate with confidence
#' inverval)
#' @param events.lab heading for events table
#' @param events.digits when estimates are shown in events table (see
#' \code{wh.events}), number of digits past the decimal to show
#' @param events.lines logical; draw lines next to groups in events table
#' @param events.col logical or a vector with colors for events table text;
#' if \code{TRUE}, \code{col.ci} will be used
#' @param main title of plot
#' @param xlab,ylab x- and y-axis labels
#' @param groups.lab labels for each line in \code{events} table
#' @param xlim,ylim x- and y-axis limits
#' @param cex.axis text size for axes labels, \code{gy_test}, abd events table
#' @param gy_test logical; if \code{TRUE} the tests of group equality will
#' be shown
#' @param test_details logical; if \code{TRUE} (default), all test details
#' (test statistic, degrees of freedom, p-value) are shown; if \code{FALSE},
#' only the p-value is shown
#' @param legend.args an optional \emph{named} list of \code{\link{legend}}
#' arguments controlling the annotations when \code{gy_test = TRUE}
#' @param split optionally split plot by unique competing risks or group;
#' one of \code{FALSE} (default, no splitting), \code{"group"}, or
#' \code{"event"}
#' @param xaxis.at,yaxis.at positions for x- and y-axis labels and ticks
#' @param xaxis.lab,yaxis.lab x- and y-axis tick labels
#' @param events.at x-coordinates to show events table (default is
#' \code{xaxis.at})
#' @param groups.order order of groups in events table
#' @param extra.margin increase left margin when groups labels in events
#' table are long (note that this will be overridden by \code{mar})
#' @param mar margins; see \code{mar} section in \code{\link{par}}
#' @param add logical; if \code{TRUE}, \code{par}s are not reset; allows for
#' multiple panels, e.g., when using \code{par(mfrow = c(1, 2))}
#' @param panel.first an expression to be evaluated after the plot axes are
#' set up but before any plotting takes place
#' @param panel.last an expression to be evaluated after plotting but before
#' returning from the function
#' @param ... additional parameters (\code{font}, \code{mfrow}, \code{bty},
#' \code{tcl}, \code{cex.lab}, \code{xaxs}, etc) passed to \code{\link{par}}
#' 
#' @seealso
#' \code{\link{cuminc2}}; \code{\link{summary.cuminc2}}
#' 
#' @examples
#' tp <- within(transplant, {
#'   futime <- futime + 1e-8
#'   age50 <- factor(+(age > 50))
#'   age_cat <- cut(age, c(0, 40, 60, Inf), c('<40', '40-60', '60+'))
#' })
#' 
#' ci1 <- cuminc2(Surv(futime, event(censored)) ~ age_cat, tp)
#' plot(ci1)
#' plot(ci1, split = 'event')
#' plot(ci1, split = 'event', events = FALSE, test_details = FALSE,
#'      legend.args = list(x = 'topright', cex = 1.5, text.col = 2,
#'                         title = "Gray\'s test p-value for"))
#' 
#' ## also plots "cuminc" objects but without extra features
#' plot(ci1$cuminc)
#' 
#' ci1 <- cuminc2(Surv(futime, event(censored) == death) ~ age50, tp)
#' plot(ci1, lty.ci = c(1,1,2,2,3,3), col.ci = 1:2)
#' 
#' ci2 <- cuminc2(Surv(futime, event(censored) == death) ~ 1, tp)
#' plot(ci2, wh.events = 'est', events.digits = 2, groups.order = c(2,1,3))
#' 
#' @export

ciplot <- function(x,
                   col.ci = seq_along(xx),
                   lty.ci = par('lty'), lwd.ci = par('lwd'),
                   
                   events = TRUE, atrisk = TRUE, events.total = TRUE,
                   wh.events = c('events', 'est', 'est.sd', 'est.ci'),
                   events.lab = NULL,
                   events.digits = 3L,
                   events.lines = TRUE, events.col = FALSE,
                   
                   main = NULL,
                   xlab = 'Time', ylab = 'Probability',
                   groups.lab = names(xx),
                   
                   xlim = NULL, ylim = NULL,
                   cex.axis = par('cex.axis'),
                   gy_test = TRUE, test_details = TRUE, legend.args = list(),
                   split = FALSE,
                   xaxis.at = pretty(xlim),
                   yaxis.at = pretty(ylim),
                   xaxis.lab = xaxis.at, yaxis.lab = yaxis.at,
                   events.at = xaxis.at,
                   groups.order = seq_along(xx),
                   
                   extra.margin = 5L, mar = NULL, add = FALSE,
                   panel.first = NULL, panel.last = NULL, ...) {
  x <- if (inherits(x, 'cuminc2'))
    x
  else if (inherits(x, 'cuminc')) {
    events <- FALSE
    list(cuminc = x, cuminc2 = NULL)
  } else assert_class(x, c('cuminc2', 'cuminc'))
  
  split <- tryCatch(match.arg(split, c('event', 'group')),
                    error = function(e) FALSE)
  
  if (!identical(split, FALSE)) {
    args <- as.list(match.call())
    args$split <- FALSE
    ## dont count censoring multiple times
    # args$atrisk <- FALSE
    
    sp <- split_cuminc(x, split)
    for (ii in sp) {
      args$x <- ii
      do.call('ciplot', args[-1L])
    }
    
    return(invisible(NULL))
  }
  
  ci <- x[['cuminc']]
  xx <- ci[names(ci) != 'Tests']
  ng <- length(xx)
  ng <- max(ng, 1L)
  cc <- col.ci
  
  groups.lab <- rep_len(groups.lab, ng)
  
  col.ci <- if (is.null(col.ci))
    seq_along(xx) else rep(col.ci, length.out = ng)
  ## must use rep instead of rep_len for names
  
  ## if cols are to be mapped by name, override col.ci
  if (!is.null(names(cc))) {
    cc <- cc[match(groups.lab, names(cc))]
    if (!anyNA(cc))
      col.ci <- col.atrisk <- cc
    else warning(
      sprintf('Mismatches in names(col.ci) - %s\n and groups.lab - %s\n',
              toString(shQuote(names(col.ci))),
              toString(shQuote(groups.lab))),
      sprintf('\nTry col.ci = c(%s)', catlist(
        setNames(as.list(shQuote(col.ci)), shQuote(groups.lab)))),
      call. = FALSE
    )
  }
  
  ## -- what is the point of this
  ## run thru tcol to convert integers and strings? idk but running a string
  ## with alpha trans removes the alpha so skip if already a hex with alpha
  # if (!any(grepl('(?i)#[a-z0-9]{8}', col.ci)))
  #   col.ci <- tcol(col.ci)
  
  lwd.ci <- rep_len(lwd.ci, ng)
  lty.ci <- rep_len(lty.ci, ng)
  
  wh.events <- match.arg(wh.events)
  
  if (!identical(events.total, FALSE)) {
    total.at <- events.total
    events.total <- TRUE
  }
  
  op <- par(no.readonly = TRUE)
  if (!add)
    on.exit(par(op))
  
  par(mar = c(4 + ng * events + atrisk,
              4 + pmax(4, extra.margin) - 3 * !events,
              2,
              2 + 3 * events * events.total))
  par(...)
  if (!is.null(mar))
    par(mar = mar)
  
  time <- unlist(sapply(xx, '[[', 'time'))
  if (is.null(xlim))
    xlim <- c(0, max(time))
  if (is.null(ylim))
    ylim <- c(0, 1)
  
  plot(
    xx[[1L]][[1L]], xx[[1L]][[2L]], type = 'n', xlim = xlim, ylim = ylim,
    ann = FALSE, axes = FALSE, panel.first = panel.first,
    panel.last = {
      box(bty = par('bty'))
      axis(1L, xaxis.at, FALSE, lwd = 0, lwd.ticks = 1)
      axis(1L, xaxis.at, xaxis.lab, FALSE, -0.5, cex.axis = cex.axis)
      axis(2L, yaxis.at, yaxis.lab, las = 1L, cex.axis = cex.axis)
      title(xlab = xlab, line = 1.5, adj = 0.5)
      title(ylab = ylab, main = main)
    })
  
  ## ci lines
  u0 <- u1 <- par('usr')
  u1[2L] <- xlim[2L] * 1.01
  do.call('clip', as.list(u1))
  for (ii in seq.int(ng))
    lines(xx[[ii]][[1L]], xx[[ii]][[2L]], ...,
          lty = lty.ci[ii], col = col.ci[ii], lwd = lwd.ci[ii])
  do.call('clip', as.list(u0))
  
  if (events) {
    usr <- par('usr')
    
    if (atrisk) {
      ng <- ng + 1L
      groups.order <- c(groups.order, ng + 1L)
      groups.lab   <- c(groups.lab, 'At-risk')[seq.int(ng)]
      
      col.ci <- c(col.ci, NA)
      lty.ci <- c(lty.ci, if (is.character(lty.ci)) 'blank' else 0)
      lwd.ci <- c(lwd.ci, 0)
    }
    
    events.at <- events.at[events.at < usr[2L]]
    
    ## set colors for lines of text
    col.events <- if (isTRUE(events.col))
      col.ci else if (!identical(events.col, FALSE) &&
                      length(events.col) == ng)
        events.col else rep_len(1L, ng)
    
    ## at-risk label
    if (!(identical(events.lab, FALSE))) {
      events.lab <- if (is.null(events.lab))
        switch(
          wh.events,
          events = 'Cumulative events',
          est    = 'Estimate',
          est.sd = 'Estimate +/- Std. dev',
          est.ci = 'Estimate [LCI, UCI]'
        )
      else events.lab
      
      mtext(events.lab, side = 1L, at = usr[1L], line = 1.5,
            adj = 1, col = 1L, las = 1L, cex = cex.axis, font = 2L)
    }
    
    ## labels for each row in at-risk table
    group.name.pos <- diff(usr[1:2]) / -8
    padding  <- abs(group.name.pos / 8)
    line.pos <- seq.int(ng)[order(groups.order)] + 2L
    
    mtext(groups.lab, side = 1L, line = line.pos, adj = 1, las = 1L,
          col = col.events, at = group.name.pos, cex = cex.axis, font = 2L)
    
    ## labels for total events
    if (events.total) {
      at <- if (isTRUE(total.at))
        abs(group.name.pos) + usr[2L] else total.at
      total.lab <- 'Total events'
      n.events <- x$cuminc2
      n.events <- n.events[n.events$status %ni% n.events$cencode, ]
      n.events <- with(droplevels(n.events), table(group, status))
      if (atrisk)
        n.events <- c(n.events, NA)
      mtext(c(n.events), side = 1L, at = at, line = line.pos, adj = 1,
            col = col.events, las = 1L, cex = cex.axis, font = 2L)
      mtext(total.lab, side = 1L, at = at, line = 1.5, adj = 1,
            col = 1L, las = 1L, cex = cex.axis, font = 2L)
    }
    
    ## draw matching lines for n events
    axpad <- 4
    if (events.lines)
      for (ii in seq.int(ng))
        ## mess with axpad here to adjust the length of the events.line
        axis(1L, c(group.name.pos + padding, 0 - axpad * padding), xpd = NA,
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
    
    if (atrisk)
      d2 <- c(d2, list('At-risk' = data.frame(
        time = as.numeric(colnames(ss)),
        n.risk = summary(x, times = events.at)$total_atrisk,
        strata = 'At-risk',
        stringsAsFactors = FALSE
      )))
    
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
            font = if (atrisk & ii == ng) 2L else 1L,
            las = 1L, line = line.pos[ii], adj = 1, col = col.events[ii])
    }
  }
  
  ## add gray text in upper right corner
  if (gy_test) {
    txt <- tryCatch(gy_text(ci, details = test_details),
                    error = function(e) NULL)
    if (is.null(txt))
      NULL
    else {
      largs <- list(
        x = 'topleft', legend = paste0(names(txt), ': ', txt),
        bty = 'n', cex = cex.axis
      )
      if (!islist(legend.args))
        legend.args <- list()
      do.call('legend', modifyList(largs, legend.args))
    }
  }
  
  panel.last
  
  invisible(NULL)
}

#' @rdname ciplot
#' @export
plot.cuminc2 <- ciplot

#' ciplot_by
#' 
#' This function helps create stratified \code{\link{ciplot}}s quickly with
#' panel labels and Gray's tests for each plot.
#' 
#' @param rhs the right-hand side of the formula
#' @param event character string indicating the event; see details
#' @param data data frame to use
#' @param by optional character string of stratification variable
#' @param single logical; if TRUE, each level of by will be drawn in a separate window
#' @param cencode unique value of \code{event} denoting censoring
#' @param gy_test one of \code{FALSE} (no test performed), a numeric value
#' (passed to \code{\link[cmprsk]{cuminc}} as the \code{rho} value), or
#' \code{TRUE} (default, \code{rho = 0})
#' @param main title(s) of plot(s)
#' @param ylab y-axis label
#' @param sub sub-title displayed in upper left corner; should be a character
#' vector with length equal to the number of panels (i.e., the number of
#' unique values of \code{by} or length one if \code{by} was not given)
#' @param groups_lab events table group labels; should be a character vector
#' with length equal to the number of groups
#' @param fig_lab figure panel labels; should be a character vector with
#' length equal to the number of panels (i.e., the number of unique values of
#' \code{by} or length one if \code{by} was not given)
#' @param col.ci color for individual curves or for all curves in a plot if
#' \code{by} is given and \code{map.col = TRUE}; if \code{col.ci} is a named
#' vector which matches the group labels, then colors are mapped to the
#' corresponding groups; see \code{\link{ciplot}}
#' @param map.col logical; if \code{TRUE}, \code{col.ci} will be the color
#' of all curves in each plot (only used when \code{by} is non-missing)
#' @param time character string of the time variable (optional)
#' @param add logical; if \code{FALSE} (default), resets graphical parameters
#' to settings before \code{ciplot_by} was called; set to \code{TRUE} for
#' adding to existing plots
#' @param plot logical; if \code{FALSE}, no plot is created but a list with
#' \code{\link{cuminc2}} object(s) is returned
#' @param ... additional arguments passed to \code{\link{ciplot}} or graphical
#' parameters subsequently passed to \code{\link{par}}
#' 
#' @return
#' Invisibly returns a list of \code{\link{cuminc2}} object(s) used to
#' generate plot(s). If \code{by} was used, there will be a list element for
#' each unique value.
#' 
#' @seealso
#' \code{\link{ciplot}}; \code{\link{cuminc2}}; \code{\link[cmprsk]{cuminc}}
#' 
#' @examples
#' ## basic usage: Surv(time, event) ~ 1
#' ciplot_by(time = 'futime', event = 'event', data = transplant)
#' 
#' ## with groups: Surv(time, event) ~ group
#' ciplot_by('sex', time = 'futime', event = 'event',
#'           data = transplant)
#' 
#' ciplot_by('sex', time = 'futime', event = 'event',
#'           data = transplant, by = 'abo', single = FALSE)
#' 
#' par(mfrow = c(1, 2))
#' ciplot_by('sex', time = 'futime', event = 'event',
#'           data = transplant, by = 'sex', xlim = c(0, 1500),
#'           single = FALSE, events.total = 2100)
#' 
#' @export

ciplot_by <- function(rhs = '1', event, data, by = NULL, single = TRUE,
                      cencode = NULL, gy_test = TRUE,
                      main = NULL, ylab = NULL, sub = NULL, groups_lab = NULL,
                      fig_lab = NULL, col.ci = NULL, map.col = FALSE,
                      time = NULL, add = FALSE, plot = TRUE, ...) {
  if (is.logical(gy_test)) {
    rho <- 0
  } else {
    rho <- if (is.numeric(gy_test))
      gy_test else 0
    gy_test <- TRUE
  }
  
  cencode <- cencode %||%
    grep('(?i)0|censor', data[, event], value = TRUE)[1L]
  cencode <- as.character(cencode)
  if (is.na(cencode))
    stop(
      '\'cencode\' not guessed correctly in ', event,
      ' variable\n  options are ', toString(shQuote(unique(data[, event]))),
      '\n  use \'cencode\' argument to specify censoring indicator'
    )
  
  form  <- if (!is.null(time))
    sprintf('Surv(%s, %s(%s)) ~ %s', time, event, cencode, rhs) else
      sprintf('Surv(%s_time, %s_ind(%s)) ~ %s', event, event, cencode, rhs)
  form  <- as.formula(form)
  
  op <- par(no.readonly = TRUE)
  if (plot & (!add | !single))
    on.exit(par(op))
  
  if (!is.null(by)) {
    data[, by] <- as.factor(data[, by])
    if (single) {
      add <- FALSE
      par(mfrow = c(1L, 1L))
    } else {
      add <- TRUE
      if (all(par('mfrow') == c(1L, 1L)))
        par(mfrow = n2mfrow(length(unique(data[, by]))))
    }
    sp <- split(data, droplevels(data[, by]))
    ## list names will be main title(s)
    names(sp) <- rep_len(main %||% paste(by, names(sp), sep = '='),
                         length(sp))
  } else {
    if (missing(add)) {
      add <- FALSE
      par(mfrow = c(1L, 1L))
    }
    map.col <- FALSE
    sp <- list(data)
  }
  
  ## define these before passing to loop
  mlabs <- is.null(groups_lab)
  msub  <- is.null(sub)
  fig   <- if (length(sp) > 1L & is.null(fig_lab))
    LETTERS[seq_along(sp)] else if (is.null(fig_lab)) '' else fig_lab
  fig   <- rep_len(fig, length(sp))
  ylab  <- ylab %||% sprintf('%s probability', toupper(event))
  
  if (!is.null(by) & is.null(names(col.ci))) {
    x  <- cuminc2(form, data)
    ci <- x[['cuminc']]
    xx <- ci[names(ci) != 'Tests']
    
    col.ci <- setNames(
      col.ci %||% seq_along(xx),
      if (by == rhs || is.null(groups_lab) || identical(groups_lab, FALSE))
        NULL else names(sl$strata)
    )
    if (length(sp) != length(col.ci))
      col.ci <- if (map.col)
        seq_along(sp) else rep(col.ci, length(sp))
  }
  
  sl <- lapply(seq_along(sp), function(x) {
    tryCatch({
      eval(substitute(
        cuminc2(form, data = sp[[x]], rho = rho),
        list(form = form))
      )
    }, error = function(e) {
      e$message
    })
  })
  
  l <- lapply(seq_along(sp), function(x) {
    s <- s0 <- sl[[x]]
    s$.data <- s0$.data <- sp[[x]]
    
    if (rhs == '1')
      rhs <- ''
    
    if (!plot)
      return(s0)
    
    ciplot(
      s, add = add, main = names(sp)[x], ylab = ylab, gy_test = gy_test, ...,
      col.ci = if (map.col) unname(col.ci)[x] else col.ci,
      panel.first = {
        ## sub label - top left margin (default: strata var)
        mtxt <- if (!msub)
          rep_len(sub, length(sp))[x] else rhs
        mtext(mtxt, line = 0.25, at = 0, adj = 0, font = 3L)
        
        ## figure label - top left outer margin (eg, A, B, C)
        mtext(fig[x], line = 0.25, at = 0 - par('usr')[2L] * .05,
              font = 2L, cex = 1.2)
      }
    )
    s0
  })
  names(l) <- names(sp) %||% rhs
  
  invisible(l)
}
