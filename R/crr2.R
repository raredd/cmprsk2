### formula method for crr
# crr2, print.crr2, summary.crr2, finegray2
# 
# unexported:
# dropcrr2, insert, tidy
###


#' Competing risks regression
#'
#' Regression modeling of subdistribution functions in competing risks.
#'
#' @param formula a \code{\link[=Surv]{survival object}} formula,
#'   \code{Surv(time, status(censor) == failure) ~ response}, where
#'   \code{censor} and \code{failure} are unique values of the \code{status}
#'   variable indicating the censoring and failure codes of interest; note all
#'   other unique values of \code{status} will be treated as competing risks
#' @param data a data frame in which to interpret the variables named in
#'   \code{formula}
#' @param which optional character string giving the desired outcome of
#'   interest (i.e., one of the unique values of the \code{status} variable);
#'   if given, no other \code{crr} models are returned
#' @param cox logical; if \code{TRUE}, a \code{\link{coxph}} model is fit
#'   using the event of interest with all other events treated as censored;
#' 
#'   alternatively, a formula to be passed to \code{\link{coxph}}, typically
#'   with the same RHS as \code{formula}, since this may be more desirable than
#'   setting \code{cox = TRUE}; note that this model may be fit to a different
#'   set of data depending on the missingness in variables of both models
#' @param variance logical; if \code{FALSE}, suppresses computation of the
#'   variance estimate and residuals
#' @param cengroup,failcode,cencode additional arguments passed to
#'   \code{\link[cmprsk]{crr}}; these will be guessed from \code{formula} but
#'   may be overridden
#' @param gtol,maxiter,init (optional) additional arguments passed to
#'   \code{\link[cmprsk]{crr}}
#' 
#' @return
#' A list of \code{\link{crr}} objects with some additional components:
#' 
#' \item{\code{coxph}}{a \code{\link{coxph}} model if \code{cox = TRUE}}
#' \item{\code{nuftime}}{a vector with the number of times each unique
#'   failure time, \code{$uftime}, was seen}
#' \item{\code{attr(, "model.frame")}}{the \code{\link{model.frame}}, i.e.,
#' \code{cov1}, used in the call to \code{\link{crr}}}
#' 
#' @seealso
#' \code{\link[cmprsk]{crr}}; \code{\link{summary.crr2}};
#' \code{\link{print.crr2}}; \code{\link{finegray2}};
#' \code{\link[survival]{finegray}}; \code{\link[survival]{coxph}}
#'
#' @examples
#' ## 'formula' as one would pass to survival::coxph with additional
#' ## indications for the censoring and failure codes
#' 
#' tp <- within(transplant, {
#'   event_ind <- as.integer(event) - 1L
#'   year <- NULL
#' })
#' 
#' ## these are equivalent ways to call 'crr2'
#' crr2(Surv(futime, event_ind(0) == 1) ~ age + sex + abo, tp)
#' crr2(Surv(futime, event_ind(0) %in% 1) ~ age + sex + abo, tp)
#' crr2(Surv(futime, event(censored) == death) ~ age + sex + abo, tp)
#' 
#' 
#' ## variables can be created or altered as usual
#' crr2(Surv(futime, event(censored) %in% death) ~ cut(age, 3), tp)
#' crr2(Surv(futime, event(censored) %in% death) ~ age + I(age ^ 2), tp)
#' crr2(Surv(futime, event(censored) %in% death) ~ relevel(abo, 'O'), tp)
#' 
#' form <- Surv(futime, event(censored) == death) ~ age + sex + abo
#' crr2(form, tp[tp$age > 55, ], cox = TRUE, variance = FALSE, which = 'death')
#' 
#' 
#' ## use the summary method to compare models easily
#' crrs <- crr2(form, data = transplant,
#'              cox = Surv(futime, event == 'death') ~ age + sex + abo)
#' summary(crrs)
#' 
#' library('htmlTable')
#' summary(crrs, html = TRUE, n = TRUE)
#' 
#' 
#' ## use the call from a crr2 object to run cmprsk::crr
#' cl <- crr2(Surv(futime, event(censored) == death) ~ age + sex + abo,
#'            na.omit(transplant))
#' cl[[1]]
#' eval(cl[[1]]$call)
#'
#'
#' ## example from cmprsk::crr
#' 
#' # simulated data to test
#' set.seed(10)
#' ftime   <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(600), 200, dimnames = list(NULL, c('x1', 'x2', 'x3')))
#' dat <- data.frame(ftime = ftime, fstatus = fstatus, cov)
#' 
#' ## cmprsk::crr
#' (z1 <- crr(ftime, fstatus, cov, failcode = 1, cencode = 0))
#' 
#' ## cmprsk2::crr2
#' (z2 <- crr2(Surv(ftime, fstatus(0) == 1) ~ ., dat))
#' 
#' summary(z1)
#' summary(z2)
#' summary(z2[[1]])
#' 
#' z1$call <- z2[[1]]$call <- NULL
#' all.equal(z1[names(z1)], z2[[1]][names(z1)])
#' # [1] TRUE
#'
#' @export

crr2 <- function(formula, data, which = NULL, cox = FALSE, variance = TRUE,
                 cengroup = NULL, failcode = NULL, cencode = NULL,
                 gtol = 1e-06, maxiter = 10, init = NULL) {
  if (!is.crr2(formula))
    stop(
      'Invalid formula - use the following structure or see ?crr2:\n',
      '\tSurv(time, status(censor) == failure) ~ response',
      call. = FALSE
    )
  
  form <- parse_formula(formula)
  name <- substitute(data)
  Name <- if (length(name) > 1L)
    as.list(name)[[2L]] else name
  
  lhs <- form$lhs
  rhs <- form$rhs
  
  # data.frame class uses "drop = TRUE" by default in single-bracket subsetting
  # but tbl_df class uses "drop = FALSE" by default. 
  # Added "drop = TRUE" argument throughout code since many of these single-bracket 
  # subsets are expecting a vector; will make code compatible with tbl_df. 
  status   <- levels(as.factor(data[, form$fstatus, drop = TRUE]))
  # cencode  <- cencode  %||% form[[2L]][1L]
  # failcode <- failcode %||% form[[2L]][2L]
  cencode  <- cencode  %||% form$cencode
  failcode <- failcode %||% form$failcode
  
  if (length(wh <- setdiff(c(lhs, rhs), names(data))))
    stop(
      sprintf('Variable%s %s not found in %s',
              ifelse(length(wh) > 1L, 's', ''),
              toString(shQuote(wh)), shQuote(name)),
      call. = FALSE
    )
  if (any(wh <- c(cencode, failcode) %ni% data[, form$fstatus, drop = TRUE]))
    warning(
      sprintf('%s not found in %s[, %s]',
              toString(shQuote(c(cencode, failcode)[wh])),
              deparse(name), shQuote(form$fstatus)),
      call. = FALSE
    )
  if (length(wh <- unique(data[, form$fstatus, drop = TRUE])) <= 2L)
    warning(
      gsub('\\n\\s{2,}', ' ',
           sprintf('Only %s found in in %s[, %s]\n Typically, censoring,
                   event of interest, and >=1 competing event are used',
                   toString(shQuote(wh)), deparse(name), shQuote(form$fstatus))
      ),
      call. = FALSE
    )
  
  ## all events of interest minus censored
  crisks <- if (!is.null(which))
    which else c(failcode, setdiff(status, c(cencode, failcode)))
  stopifnot(
    length(crisks) >= 1L
  )
  
  cengroup <- cengroup %||% rep_len(1L, nrow(data))
  
  odata <- data
  idx <- !complete.cases(data[, c(lhs, rhs)])
  if (any(!!(n <- as.integer(sum(idx))))) {
    message(n, ' observations removed due to missingness', domain = NA)
    cengroup <- cengroup[!idx]
    data <- data[!idx, ]
    assign(as.character(Name), data)
  }
  
  formula <- sprintf('Surv(%s, %s == %s) ~ %s',
                     lhs[1L], lhs[2L], shQuote(failcode),
                     paste(deparse(formula[[3L]]), collapse = ''))
  formula <- as.formula(formula)
  
  ## add model.frame for use in other methods
  mf <- model.frame(formula, data)[, -1L, drop = FALSE]
  mm <- model.matrix(formula, data)[, -1L, drop = FALSE]
  init <- rep_len(if (is.null(init)) 0L else init, ncol(mm))
  
  crrs <- lapply(crisks, function(x) {
    ftime   <- form$ftime
    fstatus <- form$fstatus
    
    ## substitute to get a more helpful call -- can run crr directly
    call <- substitute(
      crr(data[, ftime, drop = TRUE], data[, fstatus, drop = TRUE],
          cov1 = model.matrix(formula, data)[, -1L, drop = FALSE],
          cencode = cencode, failcode = x, variance = variance,
          cengroup = cengroup, gtol = gtol, maxiter = maxiter, init = init),
      list(data = name, ftime = ftime, fstatus = fstatus,
           formula = call('~', formula[[3L]]), cencode = cencode,
           x = x, variance = variance, cengroup = cengroup,
           gtol = gtol, maxiter = maxiter, init = init)
    )
    call$cengroup <- if (length(un <- unique(cengroup)) == 1L)
      substitute(rep(un, nrow(data)), list(un = un, data = name))
    else call$cengroup
    
    fit <- substitute(
      crr(data[, ftime, drop = TRUE], data[, fstatus, drop = TRUE], cov1 = mm,
          cencode = cencode, failcode = x, variance = variance,
          cengroup = cengroup, gtol = gtol, maxiter = maxiter, init = init),
      list(ftime = ftime, fstatus = fstatus, cencode = cencode,
           x = x, variance = variance, cengroup = cengroup,
           gtol = gtol, maxiter = maxiter, init = init)
    )
    fit <- eval(fit)
    fit$call <- call
    
    fit$nuftime   <- c(table(data[data[, fstatus, drop = TRUE] %in% x, ftime]))
    fit$n.missing <- n
    
    ## get n for reference group and events per model term
    fc <- +(data[, fstatus, drop = TRUE] %in% x)
    ns <- lapply(seq_along(mf), function(ii) {
      x <- mf[, ii]
      x <- if (is.factor(x) || is.character(x)) {
        x <- cbind(table(x), table(x[fc == 1L]))
        rownames(x)[1L] <- sprintf('Reference: %s', rownames(x)[1L])
        x
      } else {
        x <- cbind(length(x), sum(fc))
        rownames(x) <- colnames(mf)[ii]
        x
      }
    })
    names(ns) <- colnames(mf)
    
    ## add model.frame, model.matrix for use in other methods
    fit <- structure(
      fit, model.frame = mf, model.matrix = mm, failcode = fc, ns = ns,
      has_reference = sapply(ns, function(x)
        any(grepl('^Reference', rownames(x))))
    )
    
    structure(fit, class = c('crr2', 'crr'))
  })
  
  crrs <- structure(
    crrs, .Names = paste('CRR:', crisks),
    class = c('crr2', 'crr', 'crr2_list'),
    ## can probably get rid of this?
    ## add model.frame, model.matrix for use in other methods
    model.frame = mf, model.matrix = mm
  )
  
  if (identical(cox, FALSE))
    return(crrs)
  
  ## coxph model with any competing event vs censored
  formula <- if (isTRUE(cox)) {
    sprintf('Surv(%s, %s %%in%% c(%s)) ~ %s',
            form$ftime, form$fstatus, toString(shQuote(crisks)),
            paste(deparse(formula[[3L]]), collapse = ''))
  } else cox
  formula <- as.formula(formula)
  
  cph <- substitute(coxph(formula, data), list(formula = formula))
  cph <- eval(cph, data)
  cph$call$data <- name
  
  fc <- unname(cph$y[, 2L])
  mf <- model.frame(cph, data = data)[, -1L, drop = FALSE]
  ns <- lapply(seq_along(mf), function(ii) {
    x <- mf[, ii]
    x <- if (is.factor(x) || is.character(x)) {
      x <- cbind(table(x), table(x[fc == 1L]))
      rownames(x)[1L] <- sprintf('Reference: %s', rownames(x)[1L])
      x
    } else {
      x <- cbind(length(x), sum(fc))
      rownames(x) <- colnames(mf)[ii]
      x
    }
  })
  names(ns) <- colnames(mf)
  
  cph <- structure(
    cph, failcode = fc, model.frame = mf,
    model.matrix = model.matrix(cph), ns = ns,
    has_reference = sapply(ns, function(x)
      any(grepl('^Reference', rownames(x))))
  )
  
  structure(
    c(list('Cox PH' = cph), crrs),
    class = c('crr2', 'crr', 'crr2_list')
  )
}

#' \code{crr2} print method
#' 
#' @param x an object of class \code{\link{crr2}}
#' @param ... ignored
#' 
#' @seealso
#' \code{\link{crr2}}; \code{\link{summary.crr2}}
#' 
#' @export

print.crr2 <- function(x, ...) {
  if (any(class(x) %in% 'crr2_list')) {
    print(unclass(x)[seq_along(x)])
  } else print(dropcrr2(x))
  
  invisible(x)
}

#' \code{crr2} summary method
#'
#' @param object an object of class \code{\link{crr2}}
#' @param conf.int the level for a two-sided confidence interval on the
#'   coeficients; default is 0.95
#' @param n logical; if \code{TRUE}, the sample size and number of events
#'   for each variable are added to the summary
#' @param ref logical; if \code{TRUE}, rows with reference categories are
#'   added to results
#' @param html logical; if \code{TRUE}, an \code{\link{htmlTable}} will be
#'   returned; otherwise, a matrix
#' @param combine_ci logical; if \code{FALSE}, upper and lower confidence
#'   limits will be returned as separate columns; otherwise, an interval string
#'   will be returned
#' @param digits number of digits past the decimal point to keep
#' @param format_p logical; if \code{TRUE}, p-values will be formatted;
#'   otherwise, p-values will only be rounded
#' @param color_p logical; if \code{TRUE}, p-values will be formatted and
#'   colored based on significance level; see \code{cmprsk2:::color_pval}
#' @param format_n logical; if \code{TRUE}, for \code{html = TRUE} the total
#'   n is added for each total/events column and percents of total and events
#'   are shown in each row
#' @param htmlArgs for \code{html = TRUE}, a \emph{named} list of arguments
#'   passed to \code{\link[htmlTable]{htmlTable}} for additional formatting or
#'   to override defaults
#' @param ... ignored
#' 
#' @seealso
#' \code{\link{crr2}}; \code{\link{print.crr2}}
#' 
#' @examples
#' crrs <- crr2(Surv(futime, event(censored) == death) ~ age + sex + abo,
#'              data = transplant)
#' summary(crrs)
#' 
#' summary(crrs, conf.int = 0.9, digits = 3L,
#'         combine_ci = TRUE, format_p = TRUE)
#' 
#' library('htmlTable')
#' summary(crrs, html = TRUE, combine_ci = TRUE, n = TRUE)
#' 
#' summary(
#'   crrs, html = TRUE, combine_ci = TRUE, n = TRUE, ref = TRUE,
#'   htmlArgs = list(
#'     caption = 'CRR models.', rgroup = c('Age', 'Sex', 'Blood type'),
#'     rnames = c('+1 year change', 'Female', 'B', 'AB', 'O')
#'   )
#' )
#'
#' @export

summary.crr2 <- function(object, conf.int = 0.95, n = FALSE, ref = FALSE,
                         html = FALSE, combine_ci = html, digits = 2L,
                         format_p = html, color_p = html, format_n = n,
                         htmlArgs = list(), ...) {
  if (!any(class(object) %in% 'crr2_list'))
    return(summary(dropcrr2(object)))
  
  assert_class(object, 'crr2')
  oo <- object
  
  if (html && !is.loaded('htmlTable')) {
    message('The \'htmlTable\' package is not loaded', domain = NA)
    html <- FALSE
  }
  
  rFUN <- if (!html & !format_p & !combine_ci)
    round else roundr
  
  ns <- lapply(object, function(x) attr(x, 'ns'))
  nn <- nt <- lapply(ns, function(x) do.call('rbind', x))
  pp <- lapply(ns, function(x) {
    x <- lapply(x, function(y) {
      x <- prop.table(y, margin = 2L)
      round(x * 100)
    })
    do.call('rbind', x)
  })
  
  if (!ref) {
    nn <- lapply(nn, function(x) x[!grepl('^Reference', rownames(x)), ])
    pp <- lapply(pp, function(x) x[!grepl('^Reference', rownames(x)), ])
  }
  
  object <- lapply(object, tidy, conf.int = conf.int, ref = ref)
  object <- lapply(seq_along(object), function(ii) {
    o <- object[[ii]]
    o[, -ncol(o)] <- lapply(o[, -ncol(o)], function(x)
      rFUN(x, digits))
    if (n) {
      o <- `rownames<-`(
        cbind(`colnames<-`(nn[[ii]], c('Total', 'Events')), o),
        rownames(o)
      )
      if (html && format_n) {
        o <- within(o, {
          Total <-
            sprintf('%s (%s)', format(Total, big.mark = ','), pp[[ii]][, 1L])
          Events <-
            sprintf('%s (%s)', format(Events, big.mark = ','), pp[[ii]][, 2L])
        })
      }
    }
    o <- within(o, {
      p <- if (format_p) {
        if (color_p)
          color_pval(p)
        else pvalr(p, html = html)
      } else rFUN(p, digits)
      `HR (% CI)` <- sprintf('%s (%s, %s)', HR, LCI, UCI)
    })
    o[grepl('^Reference', rownames(o)), ncol(o)] <- NA
    
    nt <- if (html && n && ii > 1L)
      'Events' else if (n) c('Total', 'Events') else NULL
    ci <- if (combine_ci)
      'HR (% CI)' else c('HR', 'LCI', 'UCI')
    
    o[, c(nt, ci, 'p')]
  })
  names(object) <- names(oo)
  
  o <- do.call('cbind', object)
  o <- as.matrix(o)
  
  ## add in ci level
  colnames(o) <- gsub('(?=%)|(?<=[UL])CI', round(conf.int * 100),
                      colnames(o), perl = TRUE)
  ## remove crr2 list labels
  colnames(o) <- gsub('.*\\.', '', colnames(o))
  
  if (html) {
    o[grepl('^NA$|>NA<', o)] <- ''
    
    if (format_n) {
      ii <- grep('^(Total|Events)', colnames(o))
      ni <- sapply(ns, function(x) sapply(x, colSums))
      ni <- c(ni[1L, 1L], ni[2L, ])
      colnames(o)[ii] <-
        sprintf('%s<br /><font size=1>n = %s (%s)</font>',
                colnames(o)[ii], format(ni, big.mark = ','),
                round(ni / ni[1L] * 100))
      colnames(o)[1L] <- gsub('\\(\\d+\\)', '(%)', colnames(o)[1L])
    }
    
    largs <- list(
      x = o, rgroup = names(ns[[1L]]),
      n.rgroup = sapply(ns[[1L]], nrow) - (!ref) * attr(oo[[1L]], 'has_reference'),
      cgroup = c(if (n) '' else NULL, names(object)),
      n.cgroup = if (n)
        c(1L, lengths(object) - (seq_along(object) == 1L)) else lengths(object),
      css.cell = 'padding: 0px 5px 0px; white-space: nowrap;'
    )
    
    structure(
      do.call('htmlTable', modifyList(largs, htmlArgs)), class = 'htmlTable'
    )
  } else {
    div <- ifelse(combine_ci, 2L, 4L) + n * 2
    o <- lapply(seq.int(ncol(o) / div), function(x) {
      x <- x * div
      o[, (x - div + 1L):x]
    })
    setNames(o, names(oo))
  }
}

#' \code{finegray2}
#' 
#' Fit multiple \code{\link[survival]{finegray}} models.
#' 
#' @param formula a standard model formula, with survival on the left and
#'   covariates on the right
#' @param data a data frame containing the variables in \code{formula}
#' @param cencode optional censor code; if missing, the censor code will be
#'   guessed from the status variable
#' @param ... additional arguments passed to \code{\link{finegray}} or
#'   \code{\link{coxph}}
#' 
#' @examples
#' ## example from survival::finegray
#' mg <- within(mgus2, {
#'   etime <- ifelse(pstat == 0, futime, ptime)
#'   event <- ifelse(pstat == 0, 2 * death, 1)
#'   event <- factor(event, 0:2, c('censor', 'pcm', 'death'))
#' })
#' 
#' pdata <- finegray(Surv(etime, event) ~ age + sex, data = mg)
#' coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex,
#'       weight = fgwt, data = pdata)
#' 
#' fg <- finegray2(Surv(etime, event) ~ age + sex, data = mg)
#' fg$pcm
#' 
#' ## use the call and data to evaluate the equivalent coxph model
#' eval(fg$pcm$call, list(fgdata = fg$pcm$fgdata))
#' 
#' @export

finegray2 <- function(formula, data, cencode, ...) {
  localFinegray <- function(..., weights, init, control, ties, singular.ok,
                            robust, model, x, y, tt, method) {
    finegray(...)
  }
  localCoxph <- function(..., etype, prefix, count, id) {
    ## extra steps to get proper call in coxph(...)$call
    m <- match.call()
    m$etype <- m$prefix <- m$count <- m$id <- NULL
    m[[1L]] <- quote(coxph)
    fit <- coxph(...)
    fit$call <- m
    fit
  }
  
  term <- sterms.inner(formula)
  name <- substitute(data)
  
  tvar <- term[1L]
  svar <- term[2L]
  rhs <- if ('.' %in% term)
    setdiff(names(data), c(svar, tvar)) else term[-(1:2)]
  
  status  <- sort(unique(data[, svar, drop = TRUE]))
  cencode <- if (missing(cencode))
    unique(grep('(?i)0|censor', status, value = TRUE)) else cencode
  crisks  <- setdiff(status, cencode)
  
  stopifnot(
    cencode %in% status,
    length(cencode) == 1L,
    length(crisks)  >= 2L
  )
  
  data[, svar] <- factor(data[, svar, drop = TRUE], c(cencode, crisks))
  
  fg <- lapply(crisks, function(x) {
    do.call(
      'localFinegray',
      list(formula = formula, data = data, etype = x, prefix = 'fg', ...)
    )
  })
  
  cform <- sprintf(
    'Surv(fgstart, fgstop, fgstatus) ~ %s', paste(rhs, collapse = ' + ')
  )
  
  fg <- lapply(fg, function(x) {
    fit <- do.call(
      'localCoxph',
      list(formula = as.formula(cform), data = x, weights = x$fgwt, ...)
    )
    fit$call$data    <- quote(fgdata)
    fit$call$weights <- quote(fgwt)
    fit$fgdata <- x
    fit
  })
  
  setNames(fg, crisks)
}

dropcrr2 <- function(x) {
  structure(x, class = setdiff(class(x), c('crr2', 'crr2_list')))
}

insert <- function(x, where = NULL, what = paste('Reference', seq_along(where))) {
  ## insert rows into data frame/matrix and fix rownames
  # cmprsk2:::insert(head(cars), c(1, 4, 6))
  if (!length(where))
    return(x)
  
  stopifnot(
    all(where <= nrow(x))
  )
  
  ii <- seq.int(nrow(x))
  oo <- sort(c(ii, where))
  wh <- which(!duplicated(oo) & ave(oo, oo, FUN = length) > 1L)
  xo <- x[oo, ]
  xo[wh, ] <- NA
  
  if (!is.null(rn <- rownames(xo))) {
    rownames(xo)[c(wh, wh + 1L)] <- c(make.unique(what), rn[wh])
    xo
  } else xo
}

tidy <- function(x, conf.int = 0.95, ref = FALSE) {
  ## clean up crr or coxph objects
  # cmprsk2:::tidy(coxph(Surv(time, status) ~ rx, colon))
  ns <- if (ref && !is.null(ns <- attr(x, 'ns'))) {
    rr <- unlist(lapply(ns, rownames))
    ii <- grep('^Reference', rr)
    ii - head(c(0L, cumsum(!!ii)), -1L)
  } else NULL
  
  if (inherits(x, 'crr2'))
    class(x) <- c('crr', class(x))
  assert_class(x, c('crr', 'coxph'))
  
  s <- summary(x, conf.int = conf.int)
  s <- data.frame(
    s$conf.int[, -2L, drop = FALSE],
    s$coef[,  -(1:4), drop = FALSE]
  )
  s <- insert(s, ns, rr[ii])
  
  setNames(s, c('HR', 'LCI', 'UCI', 'p'))
}
