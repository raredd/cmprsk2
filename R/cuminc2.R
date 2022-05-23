### formula method for cuminc
# cuminc2, print.cuminc2, summary.cuminc2, cuminc_pairs, timepoints2
# 
# unexported:
# get_events, gy_pval, gy_text, split_cuminc, pw_pval, pw_text, name_or_index
###


#' Cumulative incidence analysis
#' 
#' Estimate cumulative incidence functions from competing risks data and test
#' equality across groups.
#' 
#' @param formula a \code{\link[=Surv]{survival object}} formula,
#'   \code{Surv(time, status(censor)) ~ response}, where \code{censor} is a
#'   unique value of the \code{status} variable indicating the censoring code;
#'   note all other unique values of \code{status} will be treated as competing
#'   risks
#' @param data a data frame in which to interpret the variables named in
#'   \code{formula}
#' @param rho,cencode,subset,na.action passed to \code{\link[cmprsk]{cuminc}};
#'   the censoring indicator will be guessed from \code{formula} but may be
#'   overridden by \code{cencode}
#' 
#' @seealso
#' \code{\link{summary.cuminc2}}; \code{\link{plot.cuminc2}};
#' \code{\link[cmprsk]{cuminc}}; 
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
#' @export

cuminc2 <- function(formula, data, rho = 0, cencode = NULL,
                    subset = NULL, na.action = getOption('na.action')) {
  if (!is.cuminc2(formula))
    stop(
      'Invalid formula - use the following structure or see ?cuminc2:\n',
      '\tSurv(time, status(censor)) ~ response',
      call. = FALSE
    )
  
  name <- substitute(data)
  Name <- if (length(name) > 1L)
    as.list(name)[[2L]] else name
  na.action <- match.fun(na.action)
  
  form  <- parse_formula(formula, data, as.character(name))
  call  <- match.call()
  call$formula <- as.formula(deparse1(formula))
  
  subset <- subset %||% rep_len(TRUE, nrow(data))
  data   <- droplevels(data[subset, ])
  
  idx <- !complete.cases(data[, c(form$lhs_vars, form$rhs_vars)])
  if (any(!!(n <- as.integer(sum(idx))))) {
    message(n, ' observations removed due to missingness', domain = NA)
    data   <- data[!idx, ]
    subset <- subset[!idx]
    assign(as.character(Name), data)
  }
  nr <- nrow(data)
  
  group  <- if (is.null(form$cov1))
    rep_len(1L, nr) else
      interaction(data[, setdiff(form$rhs_vars, form$strata), drop = FALSE])
  strata <- if (is.null(form$strata))
    rep_len(1L, nr) else
      interaction(data[, form$strata, drop = TRUE])
  
  group  <- droplevels(as.factor(group))
  strata <- droplevels(as.factor(strata))
  
  ci <- cmprsk::cuminc(
    ftime = data[, form$ftime, drop = TRUE],
    fstatus = data[, form$fstatus, drop = TRUE],
    group = group, strata = strata, rho = rho,
    cencode = cencode %||% form$cencode,
    subset = rep_len(TRUE, nr), na.action = na.action
  )
  
  ci <- list(
    cuminc  = ci,
    cuminc2 = setNames(
      data.frame(
        data[, c(form$ftime, form$fstatus)], group, strata, form$cencode,
        stringsAsFactors = FALSE
      ),
      c('time', 'status', 'group', 'strata', 'cencode')
    )
  )
  attr(ci$cuminc2, 'name') <- Name
  
  structure(
    ci, call = call, class = c('cuminc2', 'cuminc')
  )
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
#' @param digits number of digits past the decimal point to keep
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
#'   cumulative incidence model(s); see \code{\link[cmprsk]{timepoints}}}
#' \item{\code{var}}{a matrix giving the corresponding variances for the
#'   cumulative incidence model(s); see \code{\link[cmprsk]{timepoints}}}
#' \item{\code{events}}{a matrix with the number of cumulative events up
#'   to and including each value in \code{times}}
#' \item{\code{total_events}}{a vector giving the total number of events of
#'   each type (excluding censored observations)}
#' \item{\code{total_groups}}{a vector giving the number by group (including
#'   censored observations)}
#' \item{\code{total_atrisk}}{a vector giving the number remaining in the
#'   risk set regardless of event or censoring}
#' 
#' @examples
#' ci <- cuminc2(Surv(futime, event(censored)) ~ sex, transplant)
#' summary(ci)
#' 
#' ci <- cuminc2(Surv(futime, event(censored)) ~ 1, transplant)
#' summary(ci, times = 0:10 * 100)
#'
#' @export

summary.cuminc2 <- function(object, times = NULL, digits = 5L, ...) {
  assert_class(object, 'cuminc2')
  
  x <- object$cuminc2
  x <- droplevels(x[!x$status %in% x$cencode, ])
  
  gr <- levels(as.factor(x$group))
  
  sp <- split(x, interaction(x$status, x$strata, drop = TRUE))
  sp <- split(x, interaction(x$status, drop = TRUE))
  
  sp <- lapply(sp, function(x) {
    x <- x[order(x$time), ]
    x[, gr] <- lapply(gr, function(y) {
      cumsum(x$group %in% y)
    })
    x
  })
  
  if (is.null(times))
    times <- head(pretty(c(0, x$time)), -1L)
  tp <- object[['cuminc']]
  tp <- timepoints(tp, times)
  tp <- lapply(tp, round, digits = digits)
  
  
  times <- colnames(tp$est)
  sp <- lapply(sp, function(x) {
    get_events(x[, gr, drop = FALSE], x$time, times)
  })
  
  res <- do.call('rbind', sp)
  rownames(res) <- paste(gr, rep(names(sp), each = length(gr)))
  colnames(res) <- times
  
  total_events <- setNames(c(table(x$group, x$status)), rownames(res))
  
  ## "atrisk"
  total_atrisk <- object$cuminc2[order(object$cuminc2$time), ]
  total_atrisk <- get_events(NULL, total_atrisk$time, times, TRUE)
  atrisk <- c(total_events) - res
  # atrisk <- rbind(atrisk, Censored = total_atrisk - colSums(atrisk))
  
  l <- list(
    events = res, total_events = total_events,
    total_groups = c(table(object$cuminc2$group)),
    atrisk = atrisk, atrisk_sum = colSums(atrisk),
    total_atrisk = total_atrisk,
    total_censored = total_atrisk - colSums(atrisk)
  )
  
  c(tp, l)
}

get_events <- function(data, time, timepoints, atrisk = FALSE) {
  timepoints <- as.numeric(timepoints)
  
  if (atrisk) {
    res <- sapply(timepoints, function(x) sum(time > x))
    return(setNames(res, timepoints))
  }
  
  res <- lapply(timepoints, function(x) {
    idx <- time <= x
    res <- if (x < 0 || !sum(idx))
      matrix(0L, 1L, ncol(data))
    else if (sum(idx))
      data[max.col(t(idx), ties.method = 'last'), ]
    else stop('get_events - check this')
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

gy_text <- function(x, ..., details = TRUE) {
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
  if (details)
    txt else setNames(pvalr(x[, 'pv'], show.p = TRUE), rownames(x))
}

split_cuminc <- function(x, wh = c('event', 'group'), ws_split = 'last') {
  ## split (ie, hack) cuminc objects by event or group levels
  c2 <- inherits(x, 'cuminc2')
  sp_str <- switch(
    ws_split,
    last  = ' (?=\\S+$)',
    first = ' (?=.+)',
    as.numeric(ws_split)
  )
  
  ci <- if (c2)
    x[['cuminc']] else x
  xx <- ci[names(ci) != 'Tests']
  gy <- !is.null(ci$Tests)
  
  mat <- do.call(
    'rbind',
    if (is.numeric(sp_str))
      sapply(names(xx), function(x) nth(x, '\\s+', sp_str))
    else strsplit(names(xx), sp_str, perl = TRUE)
  )
  gr <- unique(mat[, 1L])
  ev <- unique(mat[, 2L])
  
  wh <- match.arg(wh)
  
  switch(
    wh,
    event = {
      
      xx_ev <- lapply(ev, function(ii) {
        tmp <- list(cuminc = xx[grep(ii, names(xx), fixed = TRUE)])
        
        if (gy)
          tmp$cuminc$Tests <- ci$Tests[grep(ii, rownames(ci$Tests)), , drop = FALSE]
        
        if (c2) {
          cc  <- as.character(x$cuminc2$cencode[1L])
          dat <- x$cuminc2[x$cuminc2[, 'status'] %in% c(ii, cc), ]
          tmp <- c(tmp, list(cuminc2 = droplevels(dat)))
        }
        
        if (c2)
          structure(tmp, class = 'cuminc2') else tmp
      })
      
      setNames(xx_ev, ev)
    },
    
    group = {
      xx_gr <- lapply(gr, function(ii) {
        tmp <- list(cuminc = xx[grep(ii, names(xx), fixed = TRUE)])
        
        if (c2) {
          dat <- x$cuminc2[x$cuminc2[, 'group'] %in% ii, ]
          tmp <- c(tmp, list(cuminc2 = droplevels(dat)))
        }
        
        if (c2)
          structure(tmp, class = 'cuminc2') else tmp
      })
      
      setNames(xx_gr, gr)
    }
  )
}

#' Pairwise \code{cuminc} comparisons
#' 
#' Evaluate pairwise group differences in cumulative incidence with
#' \code{\link[cmprsk]{cuminc}}.
#' 
#' @param object a \code{\link{cuminc2}} formula or object of class
#'   \code{\link{cuminc2}}
#' @param data a data frame to use (ignored if \code{formula}) is a
#'   \code{\link{cuminc2}} object
#' @param rho,cencode passed to \code{\link[cmprsk]{cuminc}}
#' @param method p-value correction method (default is \code{'holm'};
#'   see \code{\link{p.adjust}}
#' @param digits number of digits past the decimal point to keep
#' 
#' @examples
#' ## these two are equivalent since there are only two levels
#' form <- Surv(futime, event(censored)) ~ sex
#' cuminc_pairs(form, transplant)
#' cuminc2(form, transplant)$cuminc$Tests
#' 
#' 
#' ## more useful in the case of three or more levels
#' form <- Surv(futime, event(censored)) ~ abo
#' cuminc_pairs(form, transplant)$p.value
#' cuminc2(form, transplant)$cuminc$Tests
#' 
#' 
#' tp <- transplant[transplant$abo %in% c('A', 'B', 'O'), ]
#' tp <- tp[tp$event %in% c('censored', 'ltx', 'death'), ]
#' ci <- cuminc2(form, tp)
#' ci[[1L]]$Tests
#' ciplot(ci)
#' 
#' cuminc_pairs(ci)$p.value$ltx
#' 
#' 
#' ## unexported functions useful for extracting/formatting
#' cmprsk2:::pw_pval(ci, which = 2)
#' cmprsk2:::pw_text(ci, which = 'ltx')
#' 
#' @export

cuminc_pairs <- function(object, data = NULL, rho = 0, cencode = NULL,
                         method = p.adjust.methods, digits = 3L) {
  stopifnot(
    inherits(object, c('cuminc2', 'formula'))
  )
  
  ## pairwise gray tests
  pwgray <- function(i, j, k) {
    force(k)
    data <- data[data[, 'group', drop = TRUE] %in% c(unq[i], unq[j]), ]
    form <- sprintf('Surv(time, status(%s)) ~ group', form$cencode)
    form <- as.formula(form)
    ci <- cuminc2(form, data, rho, cencode)
    
    tryCatch(
      gy_pval(ci, TRUE)[k, 'stat'],
      error = function(e) {
        if (grepl('subscript', e$message))
          NA else stop(e)
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
  unq  <- levels(as.factor(data[, 'group', drop = TRUE]))
  ngy  <- nrow(object[['cuminc']]$Tests)
  crs  <- setdiff(levels(as.factor(data[, 'status', drop = TRUE])), form$cencode)
  
  nn <- outer(as.character(unq), as.character(unq), Vectorize(function(x, y)
    nrow(data[data[, 'group', drop = TRUE] %in% c(x, y), ])))
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
  
  structure(
    list(n = nn, chi.sq = chisq, p.value = p.value),
    class = 'cuminc_pairs'
  )
}

pw_pval <- function(object, details = FALSE, data = NULL, ...,
                    method = 'none', which = NULL) {
  object <- if (inherits(object, 'cuminc_pairs'))
    object
  else if (inherits(object, 'cuminc2')) {
    cuminc_pairs(object, ..., digits = 10L)
  } else if (inherits(object, c('formula', 'call'))) {
    object <- eval(
      substitute(cuminc2(form, data = data, ...), list(form = object))
    )
    cuminc_pairs(object)
  } else stop('pw_pval - Invalid object', call. = FALSE)
  
  stopifnot(
    inherits(object, 'cuminc_pairs'),
    !is.null(which)
  )
  
  m <- object$p.value
  m <- m[[name_or_index(which, names(m))]]
  p <- m[lower.tri(m)]
  p <- p.adjust(p, method = method)
  n <- sprintf('%s vs %s', colnames(m)[col(m)[lower.tri(m, FALSE)]],
               rownames(m)[row(m)[lower.tri(m, FALSE)]])
  
  setNames(p, n)
}

pw_text <- function(formula, data, ..., details = TRUE, pFUN = NULL,
                    method = 'none', which = NULL) {
  pFUN <- if (is.null(pFUN) || isTRUE(pFUN))
    function(x) pvalr(x, show.p = TRUE)
  else if (identical(pFUN, FALSE))
    identity
  else {
    stopifnot(is.function(pFUN))
    pFUN
  }
  
  obj <- pw_pval(object = formula, data = data, method = method,
                 which = which, ...)
  
  sprintf('%s: %s', names(obj), pFUN(obj))
}

name_or_index <- function(x, y = NULL) {
  ## return integer vector where x occurs in y
  # cmprsk2:::name_or_index(c('1', '3', 'e'))
  # cmprsk2:::name_or_index(c('a', 'c', 'e'), letters)
  # table is given priority over integer, eg, idx = 27 instead of 4
  # cmprsk2:::name_or_index(c('a', '4', 'e'), c(letters, '4'))
  suppressWarnings(
    ix <- as.integer(x)
  )
  
  if (!is.null(y)) {
    iy <- match(x, y)
    replace(iy, is.na(iy), ix[is.na(iy)])
  } else ix
}

#' \code{timepoints} formatter
#' 
#' Format output from \code{\link[cmprsk]{timepoints}} into a single matrix
#' with estimate +/- standard deviation.
#' 
#' @param w,times arguments passed to \code{\link{timepoints}}
#' @param digits number of digits past the decimal point to keep
#' @param sd logical; if \code{FALSE}, the standard deviation will not be
#'   shown with the estimate
#' @param ci logical; not implemented
#' @param html logical; if \code{TRUE}, an html-friendly format is returned;
#'   the print method for \code{timepoints2} will use \code{\link{htmlTable}}
#'   if \code{html = TRUE}
#' @param htmlArgs for \code{html = TRUE}, a \emph{named} list of arguments
#'   passed to \code{\link[htmlTable]{htmlTable}} for additional formatting or
#'   to override defaults
#' @param ... ignored
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
#' gg <- factor(sample(1:3, 100, replace = TRUE), 1:3, c('a', 'b', 'c'))
#' cc <- sample(0:2, 100, replace = TRUE)
#' strt <- sample(1:2, 100, replace = TRUE)
#' xx <- cuminc(ss, cc, gg, strt)
#' 
#' timepoints(xx, times = 0:4)
#' timepoints2(xx, times = 0:4)
#' 
#' @export

timepoints2 <- function(w, times = NULL, digits = 3L, sd = FALSE, ci = FALSE,
                        html = FALSE, htmlArgs = list(), ...) {
  w <- if (inherits(w, 'cuminc2'))
    w[['cuminc']]
  else if (inherits(w, 'cuminc'))
    w else stop('\'w\' should be of class \'cuminc\' or \'cuminc2\'')
  
  times <- if (is.null(times)) {
    st <- sort(unlist(sapply(w, `[`, 'time')))
    pt <- pretty(st)
    pt[pt %inside% range(st)]
  } else sort(unique(times))
  
  tp  <- timepoints(w, times)
  res <- tp$est
  
  fmt <- if (!ci)
    sprintf('%%.%sf %s %%.%sf', digits, if (html) '&pm;' else '+/-', digits)
  else sprintf('%%.%sf [%%.%sf - %%.%sf]', digits, digits, digits)
  
  if (ci) {
    z <- qnorm(1 - 0.05 / 2)
    res[] <- mapply(function(x, y, z)
      sprintf(fmt, x, y, z),
      tp$est,
      pmax(tp$est - z * sqrt(tp$var), 0),
      tp$est + z * sqrt(tp$var))
  } else {
    res[] <- mapply(function(x, y)
      sprintf(fmt, x, y), tp$est, sqrt(tp$var))
    if (!sd)
      res <- gsub(' [+&].*$', '', res)
  }
  
  res <- gsub('NA.*$', '-', res)
  
  if (html) {
    largs <- list(
      x = res, css.cell = 'padding: 0px 5px 0px; white-space: nowrap;'
    )
    
    structure(
      do.call('htmlTable', modifyList(largs, htmlArgs)), class = 'htmlTable'
    )
  } else res
}
