### formula method for crr
# crr2, tidy_, summary.crr2, print.crr2
###


### TODO
# finegray2
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
#' 
#' library('htmlTable')
#' summary(crrs)
#' 
#' crr2(form, transplant, cox = TRUE, variance = FALSE)
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
#' # deviance(z)
#'
#' dat <- data.frame(ftime = ftime, fstatus = fstatus, cov)
#' (z2 <- crr2(Surv(ftime, fstatus(0) == 1) ~ ., dat, cox = FALSE))
#' z$call <- z2[[1]]$call <- NULL
#' stopifnot(identical(z, z2[[1]]))
#'
#' @export

crr2 <- function(formula, data, which = NULL, cox = FALSE, ...) {
  term <- terms.inner(formula)
  form <- rapply(term, trimwsq, how = 'list')
  name <- substitute(data)
  
  if ((!attr(term, 'dots') & any(form[[1L]] %ni% names(data))) ||
      length(form[[2L]]) != 2L ||
      any(form[[2L]] %ni% data[, form[[1L]][2L]]))
    stop('Invalid formula - use the following structure or see ?crr2:\n',
         '\tSurv(time, status(censor) == failure) ~ response', call. = FALSE)

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
  } else
    n <- 0L

  cov1 <- model.matrix(formula, data)[, -1L, drop = FALSE]
  crrs <- lapply(crisks, function(x) {
    ftime <- lhs[1L]
    fstatus <- lhs[2L]
    cl <- substitute(
      crr(data[, ftime], data[, fstatus], cov1 = cov1,
          cencode = cencode, failcode = x, ...),
      list(data = name, ftime = ftime, fstatus = fstatus,
           cencode = cencode, x = x)
    )
    fit <- eval(cl)
    # fit <- crr(data[, lhs[1L]], data[, lhs[2L]], cov1 = cov1,
    #            cencode = cencode, failcode = x, ...)
    fit$n.missing <- n
    fit
  })
  crrs <- structure(crrs, .Names = paste('CRR:', crisks), class = 'crr2')
  
  ## add model.frame for other use
  mf <- model.frame(reformulate(rhs, lhs[1L]), data)
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

  structure(c(cph, crrs), class = 'crr2')
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
#' @param object an object of class \code{"crr2"}
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
  
  if (!is.loaded('htmlTable')) {
    message('The \'htmlTable\' package is not loaded', domain = NA)
    html <- FALSE
  }
  
  if (html) {
    # pv <- grep('^p$', names(o))
    # o[pv] <- lapply(o[pv], function(x) {
    #   x <- gsub('<', '&lt;', x, fixed = TRUE)
    #   x <- gsub('>', '&gt;', x, fixed = TRUE)
    #   x
    # })
    o <- gsub('<', '&lt;', o, fixed = TRUE)
    o <- gsub('>', '&gt;', o, fixed = TRUE)
    htmlTable(o, ..., css.cell = 'padding: 0px 5px 0px;',
              cgroup = names(object), n.cgroup = lengths(object))
  } else o
}

#' @rdname crr2
#' @export
print.crr2 <- function(x, ...) {
  print(x[seq_along(x)])
  invisible(x)
}

# finegray2 <- function(formula, data, which = NULL, cox = FALSE, ...) {
#   # Surv(time, status(censored) == interest) ~ vars
#   # or
#   # Surv(time, status) ~ vars
#   
#   term <- cmprsk2:::terms.inner(formula)
#   data[, term[2L]] <- as.factor(data[, term[2L]])
#   
# }
