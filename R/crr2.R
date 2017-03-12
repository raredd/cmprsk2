### formula method for crr
# crr2, tidy_, summary.crr2, print.crr2, finegray2
###


#' Competing risks regression
#'
#' Regression modeling of subdistribution functions in competing risks.
#'
#' @param formula a formula object of form
#' \code{Surv(time, status(censor) == failure) ~ response} where \code{censor}
#' and \code{failure} are labels of the \code{status} variable denoting
#' the censor and failure of interest; note all other unique values of
#' \code{status} will be considered competing risks
#' @param data a data.frame in which to interpret the variables named in
#' \code{formula}
#' @param which optional character string giving the desired outcome of
#' interest; if given, no other models are returned
#' @param cox logical; if \code{TRUE}, a \code{\link{coxph}} model is also
#' fit modeling the event of interest
#' @param ... additional arguments passed to \code{\link[cmprsk]{crr}}
#' 
#' @return
#' A list of \code{\link{crr}} objects with some additional attributes:
#' 
#' \item{\code{$nuftime}}{a vector with the number of times each unique
#' failure time \code{$uftime} was seen}
#' \item{\code{attr(, "model.frame")}}{the \code{\link{model.frame}}, i.e.,
#' \code{cov1}, used in the call to \code{\link{crr}}}
#' 
#'
#' @seealso
#' \code{\link{summary.crr2}}; \code{\link[cmprsk]{crr}};
#' \code{\link[survival]{coxph}}
#'
#' @examples
#' crr2(Surv(futime, event(censored) == death) ~ age + sex + abo, transplant)
#' 
#' form <- Surv(futime, event(censored) == death) ~ age + sex + abo
#' crrs <- crr2(form, transplant, which = 'ltx', cox = TRUE)
#' 
#' summary(crrs)
#' library('htmlTable')
#' summary(crrs)
#' 
#' crr2(form, transplant, cox = TRUE, variance = FALSE)
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
#' ftime <- rexp(200)
#' fstatus <- sample(0:2,200,replace=TRUE)
#' cov <- matrix(runif(600),nrow=200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3')
#' print(z <- crr(ftime,fstatus,cov))
#' summary(z)
#' deviance(z)
#'
#' dat <- data.frame(ftime = ftime, fstatus = fstatus, cov)
#' (z2 <- crr2(Surv(ftime, fstatus(0) == 1) ~ ., dat, cox = FALSE))
#' z$call <- z2[[1]]$call <- NULL
#' all.equal(z[names(z)], z2[[1]][names(z)])
#'
#' @export

crr2 <- function(formula, data, which = NULL, cox = FALSE, ...) {
  term <- terms.inner(formula)
  form <- rapply(term, trimwsq, how = 'list')
  name <- substitute(data)
  
  # if ((!attr(term, 'dots') & any(form[[1L]] %ni% names(data))) ||
  #     length(form[[2L]]) != 2L ||
  #     any(form[[2L]] %ni% data[, form[[1L]][2L]]))
  #   stop('Invalid formula - use the following structure or see ?crr2:\n',
  #        '\tSurv(time, status(censor) == failure) ~ response', call. = FALSE)

  lhs  <- form[[1L]][1:2]
  rhs  <- if (attr(term, 'dots'))
    setdiff(names(data), lhs) else form[[1L]][-(1:2)]
  
  status   <- sort(unique(data[, lhs[2L]]))
  cencode  <- form[[2L]][1L]
  failcode <- form[[2L]][2L]

  formula <- sprintf('Surv(%s, %s == %s)  ~ %s',
                     lhs[1L], lhs[2L], shQuote(failcode),
                     paste(rhs, collapse = ' + '))
  formula <- as.formula(formula)

  ## all events of interest minus censored
  crisks  <- if (!is.null(which))
    which else c(failcode, setdiff(status, c(cencode, failcode)))
  stopifnot(length(crisks) >= 1L, failcode %in% status, cencode %in% status)

  idx <- !complete.cases(data[, c(lhs, rhs)])
  if (any(idx)) {
    n <- as.integer(sum(idx))
    message(n, ' observations removed due to missingness', domain = NA)
    data <- data[!idx, ]
    assign(as.character(name), data)
  } else
    n <- 0L

  ## add model.frame for other use
  mf <- model.frame(reformulate(rhs, lhs[1L]), data)
  # cov1 <- model.matrix(formula, data)[, -1L, drop = FALSE]
  
  crrs <- lapply(crisks, function(x) {
    ftime <- lhs[1L]
    fstatus <- lhs[2L]
    # fit <- crr(data[, lhs[1L]], data[, lhs[2L]], cov1 = cov1,
    #            cencode = cencode, failcode = x, ...)
    
    ## substitute to get a more helpful call -- can run crr directly
    fit <- substitute(
      crr(data[, ftime], data[, fstatus],
          cov1 = model.matrix(formula, data)[, -1L, drop = FALSE],
          cencode = cencode, failcode = x, ...),
      list(data = name, ftime = ftime, fstatus = fstatus,
           formula = update(formula, ~ .), cencode = cencode, x = x)
    )
    fit <- eval(fit)
    
    fit$nuftime     <- c(table(data[data[, fstatus] %in% x, ftime]))
    fit$n.missing   <- n
    attr(fit, 'model.frame') <- mf
    
    structure(fit, class = c('crr', 'crr2'))
  })
  crrs <- structure(
    crrs, .Names = paste('CRR:', crisks), class = c('crr2', 'crr')
  )
  
  ## add model.frame for other use
  attr(crrs, 'model.frame') <- mf
  # names(crrs) <- paste('CRR:', crisks)

  if (!cox)
    return(crrs)

  ## collapse all events to fit coxph (or adjust formula above)
  # data[, lhs[2]] <- +(data[, lhs[2]] %in% crisks)
  # data[, lhs[2]] <- +(data[, lhs[2]] %in% failcode)
  cph <- substitute(survival::coxph(formula, data), list(formula = formula))
  cph <- eval(cph)
  cph$call$data <- name
  cph <- list('Cox PH' = cph)

  structure(c(cph, crrs), class = c('crr2', 'crr'))
}

tidy_ <- function(x, conf.int = 0.95, ...) {
  ## clean up crr or coxph objects
  # tidy_(coxph(Surv(time, status) ~ rx, colon))
  assert_class(x, c('crr', 'coxph'))
  s <- summary(x, conf.int = conf.int, ...)
  
  setNames(data.frame(s$conf.int[, -2L, drop = FALSE],
                      s$coef[,  -(1:4), drop = FALSE]),
           c('HR', 'LCI', 'UCI', 'p'))
}

#' \code{crr2} summary method
#'
#' @param object an object of class \code{\link{crr2}}
#' @param conf.int level of confidence
#' @param ... additional arguments affecting the summary produced
#' @param html logical; if \code{TRUE}, an \code{\link{htmlTable}} will be
#' returned; otherwise, a matrix
#' @param combine_ci logical; if \code{FALSE}, upper and lower confidence
#' limits will be returned as separate columsn; otherwise, an interval string
#' will be returned
#' @param digits integer value indicating number of digits to print
#' @param format_p logical; if \code{TRUE}, p-values will be formatted;
#' otherwise, p-values will only be rounded
#'
#' @export

summary.crr2 <- function(object, conf.int = 0.95, ..., html = TRUE,
                         combine_ci = FALSE, digits = 2L, format_p = TRUE) {
  stopifnot(inherits(object, 'crr2'))
  
  object <- lapply(object, tidy_, conf.int = conf.int, ...)
  object <- lapply(object, function(o) {
    o[, -ncol(o)] <- lapply(o[, -ncol(o)], function(x)
      sprintf('%.*f', digits, x))
    within(o, {
      p <- if (format_p)
        format.pval(p, 1, .01) else sprintf('%.*f', digits, p)
      `HR (% CI)` <- sprintf('%s (%s, %s)', HR, LCI, UCI)
    })[, if (combine_ci) c('HR (% CI)', 'p') else -(ncol(o) + 1L)]
  })
  
  o <- do.call('cbind', object)
  o <- as.matrix(o)
  
  ## add in ci level
  colnames(o) <- gsub('(?=%)|(?<=[UL])CI', round(conf.int * 100),
                      colnames(o), perl = TRUE)
  ## remove crr2 list labels
  colnames(o) <- gsub('.*\\.', '', colnames(o))
  
  if (html && !is.loaded('htmlTable')) {
    message('The \'htmlTable\' package is not loaded', domain = NA)
    html <- FALSE
  }
  
  if (html) {
    o <- gsub('<', '&lt;', o, fixed = TRUE)
    o <- gsub('>', '&gt;', o, fixed = TRUE)
    
    ## bug in htmlTable v1.9 with class == c('html', ...)
    structure(
      htmlTable(o, ..., cgroup = names(object), n.cgroup = lengths(object),
                css.cell = 'padding: 0px 5px 0px; white-space: nowrap;'),
      class = 'htmlTable'
    )
  } else {
    o <- lapply(seq.int(ncol(o) / 4), function(x) {
      x <- x * 4
      o[, (x - 4 + 1):x]
    })
    setNames(o, names(object))
  }
}

#' \code{crr2} print method
#' 
#' @param x an object of class \code{\link{crr2}}
#' @param ... ignored
#' 
#' @seealso
#' \code{\link{crr2}}
#' 
#' @export
print.crr2 <- function(x, ...) {
  if (inherits(x, 'crr2')) {
    print(unclass(x)[seq_along(x)])
  } else print(x)
  invisible(x)
}

#' \code{finegray2}
#' 
#' Fit multiple \code{\link[survival]{finegray}} models.
#' 
#' @param formula a standard model formula, with survival on the left and
#' covariates on the right
#' @param data a data frame containing the variables in \code{formula}
#' @param cencode optional censor code; if missing, the censor code will be
#' guessed from the status variable
#' @param ... additional arguments passed to \code{\link{finegray}} or
#' \code{\link{coxph}}
#' 
#' @examples
#' ## example from survival::finegray
#' dd <- within(mgus2, {
#'   etime <- ifelse(pstat == 0, futime, ptime)
#'   event <- ifelse(pstat == 0, 2 * death, 1)
#'   event <- factor(event, 0:2, c('censor', 'pcm', 'death'))
#' })
#' 
#' pdata <- finegray(Surv(etime, event) ~ age + sex, data = dd)
#' coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex,
#'       weight = fgwt, data = pdata)
#' 
#' fg <- finegray2(Surv(etime, event) ~ age + sex, data = dd)
#' fg$pcm
#' 
#' eval(fg$pcm$call, list(fgdata = fg$pcm$fgdata))
#' 
#' @export

finegray2 <- function(formula, data, cencode, ...) {
  localFinegray <- function(..., weights, init, control, ties, singular.ok,
                            robust, model, x, y, tt, method)
    finegray(...)
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
  
  status  <- sort(unique(data[, svar]))
  cencode <- if (missing(cencode)) {
    unique(grep('(?i)0|censor', status, value = TRUE))
  } else {
    cencode
  }
  crisks <- setdiff(status, cencode)
  
  stopifnot(
    cencode %in% status,
    length(cencode) == 1L,
    length(crisks)  >= 2L
  )
  
  data[, svar] <- factor(data[, svar], c(cencode, crisks))
  
  fg <- lapply(crisks, function(x)
    do.call('localFinegray', list(formula = formula, data = data, etype = x,
                                  prefix = 'fg', ...)))
  
  cform <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s',
                   paste(rhs, collapse = ' + '))
  
  fg <- lapply(fg, function(x) {
    fit <- do.call('localCoxph', list(formula = as.formula(cform), data = x,
                                      weights = x$fgwt, ...))
    fit$call$data    <- quote(fgdata)
    fit$call$weights <- quote(fgwt)
    fit$fgdata <- x
    fit
  })
  
  setNames(fg, crisks)
}
