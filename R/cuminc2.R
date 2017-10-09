### formula method for cuminc
# cuminc2, ciplot, plot.cuminc2, ciplot_by, print.cuminc2, summary.cuminc2,
# cuminc_pairs
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
  
  form  <- parse_formula(formula, data)
  fname <- as.formula(deparse(formula))
  call  <- match.call()
  call$formula <- fname
  
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
  attr(ci$cuminc2, 'name') <- Name
  
  structure(
    ci,
    call = call,
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

ciplot <- function(x,
                   main = NULL,
                   xlab = 'Time', ylab = 'Probability',
                   groups.lab = names(xx),
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
                   groups.order = seq_along(xx),
                   extra.margin = 5L,
                   panel.first = NULL, panel.last = NULL, 
                   add = FALSE, mar = NULL, ...) {
  x <- if (inherits(x, 'cuminc2'))
    x
  else if (inherits(x, 'cuminc')) {
    events <- FALSE
    list(cuminc = x, cuminc2 = NULL)
  } else assert_class(x, c('cuminc2', 'cuminc'))
  
  ci <- x[['cuminc']]
  xx <- ci[names(ci) != 'Tests']
  ng <- length(xx)
  ng <- max(ng, 1L)
  cc <- col.ci
  
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
  
  ## run thru tcol to convert integers and strings? idk but
  ## running a string with alpha trans resets to no alpha/trans
  ## so skip if already a hex color with alpha/trans
  if (!any(grepl('(?i)#[a-z0-9]{8}', col.ci)))
    col.ci <- tcol(col.ci)
  
  lwd.ci <- rep_len(lwd.ci, ng)
  lty.ci <- rep_len(lty.ci, ng)
  
  wh.events <- match.arg(wh.events)
  
  op <- par(no.readonly = TRUE)
  if (!add)
    on.exit(par(op))
  
  par(mar = c(4 + ng * events,
              4 + pmax(4, extra.margin) - 3 * !events,
              2,
              2 + 3 * events))
  par(...)
  if (!is.null(mar))
    par(mar = mar)
  
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
  u0 <- u1 <- par('usr')
  u1[2L] <- xlim[2L]
  do.call('clip', as.list(u1))
  for (ii in seq.int(ng))
    lines(xx[[ii]][[1L]], xx[[ii]][[2L]], ...,
          lty = lty.ci[ii], col = col.ci[ii], lwd = lwd.ci[ii])
  do.call('clip', as.list(u0))
  
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
    line.pos <- seq.int(ng)[order(groups.order)] + 2L
    
    mtext(groups.lab, side = 1L, line = line.pos, adj = 1, las = 1L,
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
#' 
#' @examples
#' ## basic usage: Surv(time, event) ~ 1
#' ciplot_by(time = 'futime', event = 'event', data = transplant)
#' 
#' ## with groups: Surv(time, event) ~ group
#' ciplot_by('sex', time = 'futime', event = 'event',
#'           data = transplant)
#' ciplot_by('sex', time = 'futime', event = 'event',
#'           data = transplant, by = 'abo', single = FALSE)
#' 
#' @export

ciplot_by <- function(rhs = '1', event, data, by = NULL, single = TRUE,
                      cencode = NULL, gy_test = TRUE,
                      main = NULL, ylab = NULL, sub = NULL, strata_lab = NULL,
                      fig_lab = NULL, col.ci = NULL, map.col = FALSE,
                      time = NULL, add = FALSE, plot = TRUE, ...) {
  if (is.logical(gy_test)) {
    rho <- 0
  } else if (is.numeric(gy_test)) {
    rho <- gy_test
    gy_test <- TRUE
  } else {
    rho <- 0
    gy_test <- TRUE
  }
  
  cencode <- cencode %||%
    grep('(?i)0|censor', data[, event], value = TRUE)[1L]
  cencode <- as.character(cencode)
  
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
  mlabs <- is.null(strata_lab)
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
      if (by == rhs || is.null(strata_lab) || identical(strata_lab, FALSE))
        NULL else strata_names %||% names(sl$strata)
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
    
    ciplot(s, add = add, main = names(sp)[x], ylab = ylab,
           gy_test = gy_test, ...,
           col.ci = if (map.col) unname(col.ci)[x] else col.ci,
           panel.first = {
             ## sub label - top left margin (default: strata var)
             mtxt <- if (!msub)
               rep_len(sub, length(sp))[x] else rhs
             mtext(mtxt, 3, 0.25, FALSE, 0, 0, font = 3)
             
             ## figure label - top left outer margin (eg, A, B, C)
             mtext(fig[x], 3, 0.25, FALSE, 0 - par('usr')[2L] * .05,
                   font = 2, cex = 1.2)
           })
    s0
  })
  names(l) <- names(sp) %||% rhs
  
  invisible(l)
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

#' Pairwise \code{cuminc} comparisons
#' 
#' Evaluate pairwise group differences in cumulative incidence with
#' \code{\link[cmprsk]{cuminc}}.
#' 
#' @param object a \code{\link{cuminc2}} formula or object of class
#' \code{\link{cuminc2}}
#' @param data a data frame to use (ignored if \code{formula}) is a
#' \code{\link{cuminc2}} object
#' @param rho,cencode passed to \code{\link[cmprsk]{cuminc}}
#' @param method p-value correction method (default is \code{'holm'};
#' see \code{\link{p.adjust}}
#' @param digits integer indicating the number of decimal places to be used
#' 
#' @examples
#' ## these two are equivalent since there are only two levels
#' form <- Surv(futime, event(censored)) ~ sex
#' cuminc_pairs(form, transplant)
#' cuminc2(form, transplant)$cuminc$Tests
#' 
#' ## more useful in the case of three or more levels
#' form <- Surv(futime, event(censored)) ~ abo
#' cuminc_pairs(form, transplant)$p.value
#' cuminc2(form, transplant)$cuminc$Tests
#' 
#' tp <- transplant[transplant$abo %in% c('A', 'B', 'O'), ]
#' tp <- tp[tp$event %in% c('censored', 'ltx', 'death'), ]
#' ci <- cuminc2(form, tp)
#' ci[[1L]]$Tests
#' ciplot(ci)
#' 
#' cuminc_pairs(ci)$p.value$ltx
#' 
#' @export

cuminc_pairs <- function(object, data = NULL, rho = 0, cencode = NULL,
                         method = p.adjust.methods, digits = 3L) {
  pwgray <- function(i, j, k) {
    force(k)
    data <- data[data[, 'group'] %in% c(unq[i], unq[j]), ]
    form <- sprintf('Surv(time, status(%s)) ~ group', form$cencode)
    form <- as.formula(form)
    ci <- cuminc2(form, data, rho, cencode)
    
    tryCatch(gy_pval(ci, TRUE)[k, 'stat'],
             error = function(e) {
               if (grepl('subscript', e$message))
                 NA
               else stop(e)
             })
  }
  
  ## stats::pairwise.table with modifications
  pwtable <- function(compare.levels, level.names, p.adjust.method, ...) {
    ix <- setNames(seq_along(level.names), level.names)
    pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec)
      sapply(seq_along(ivec), function(k) {
        i <- ivec[k]
        j <- jvec[k]
        if (i > j)
          compare.levels(i, j, ...) else NA
      }))
    
    pp[lower.tri(pp, TRUE)] <-
      p.adjust(pp[lower.tri(pp, TRUE)], p.adjust.method)
    pp
  }
  
  if (inherits(object, 'formula')) {
    suppressMessages(
      object <- cuminc2(object, data, rho, cencode)
    )
  }
  
  method <- match.arg(method)
  
  data <- droplevels(object[['cuminc2']])
  form <- parse_formula(attr(object, 'call')$formula)
  unq  <- levels(as.factor(data[, 'group']))
  ngy  <- nrow(object[['cuminc']]$Tests)
  crs  <- setdiff(levels(as.factor(data[, 'status'])), form$cencode)
  
  nn <- outer(as.character(unq), as.character(unq), Vectorize(function(x, y)
    nrow(data[data[, 'group'] %in% c(x, y), ])))
  nn[upper.tri(nn, FALSE)] <- NA
  dimnames(nn) <- list(unq, unq)
  
  chisq <- lapply(seq.int(ngy), function(ii) {
    chi <- rbind(NA, cbind(pwtable(pwgray, unq, 'none', k = crs[ii]), NA))
    dimnames(chi) <- list(unq, unq)
    round(chi, digits)
  })
  
  p.value <- lapply(chisq, function(chi) {
    pv <- apply(chi, 1:2, function(x)
      pchisq(x, 1L, lower.tail = FALSE))
    tpv <- t(pv)
    tpv[lower.tri(pv)] <- p.adjust(pv[lower.tri(pv)], method = method)
    round(t(tpv), digits)
  })
  names(chisq) <- names(p.value) <- crs
  
  list(n = nn, chi.sq  = chisq, p.value = p.value)
}
