### cmprsk functions with improvements
# crr
###


crr <- function(ftime, fstatus, cov1 = NULL, cov2 = NULL, tf = NULL,
                cengroup = NULL, failcode = 1, cencode = 0, subset = NULL,
                na.action = na.omit, gtol = 1e-06, maxiter = 10, init = NULL,
                variance = TRUE) {
  call <- match.call()
  cov1.name <- if (is.null(cov1))
    '' else deparse1(substitute(cov1))
  cov1.vars <- cov2.vars <- NULL
  if (!is.null(cov1)) {
    cov1.vars <- colnames(as.matrix(cov1))
  }
  cov2.name <- if (is.null(cov2))
    '' else deparse1(substitute(cov2))
  if (!is.null(cov2)) {
    cov2.vars <- colnames(as.matrix(cov2))
  }
  
  d <- data.frame(
    ftime = ftime,
    fstatus = fstatus,
    cengroup = if (is.null(cengroup))
      rep(1, length(fstatus)) else cengroup
  )
  if (!is.null(cov1)) {
    cov1 <- as.matrix(cov1)
    nc1 <- ncol(cov1)
    d <- cbind(d, cov1)
  } else {
    nc1 <- 0
  }
  if (!is.null(cov2)) {
    cov2 <- as.matrix(cov2)
    nc2 <- ncol(cov2)
    d <- cbind(d, cov2)
  } else {
    nc2 <- 0
  }
  
  if (!is.null(subset)) 
    d <- d[subset, ]
  
  tmp <- nrow(d)
  d <- na.action(d)
  nmis <- 0
  
  if (nrow(d) != tmp) {
    nmis <- tmp - nrow(d)
    cat(format(nmis), 'cases omitted due to missing values\n')
  }
  
  d <- d[order(d$ftime), ]
  
  ftime <- d$ftime
  cenind <- ifelse(d$fstatus == cencode, 1, 0)
  fstatus <- ifelse(d$fstatus == failcode, 1, 2 * (1 - cenind))
  
  ucg <- sort(unique.default(d$cengroup))
  cengroup <- match(d$cengroup, ucg)
  ncg <- length(ucg)
  uuu <- matrix(0, nrow = ncg, ncol = length(ftime))
  
  for (k in 1:ncg) {
    u <- do.call(
      'survfit',
      list(formula = Surv(ftime, cenind) ~ 1,
           data = data.frame(ftime, cenind, cengroup), subset = cengroup == k)
    )
    u <- approx(
      c(min(0, u$time) - 10 * .Machine$double.eps,
        c(u$time, max(u$time) * (1 + 10 * .Machine$double.eps))),
      c(1, u$surv, 0), xout = ftime * (1 - 100 * .Machine$double.eps),
      method = 'constant', f = 0, rule = 2
    )
    uuu[k, 1:length(u$y)] <- u$y
  }
  
  uft <- sort(unique(ftime[fstatus == 1]))
  ndf <- length(uft)
  
  if (nc2 == 0) {
    cov1 <- as.matrix(d[, (1:nc1) + 3])
    np <- nc1
    npt <- 0
    cov2 <- 0
    tfs <- 0
  } else if (nc1 == 0) {
    cov2 <- as.matrix(d[, (1:nc2) + 3 + nc1])
    npt <- np <- nc2
    cov1 <- 0
    tfs <- tf(uft)
  } else {
    cov1 <- as.matrix(d[, (1:nc1) + 3])
    cov2 <- as.matrix(d[, (1:nc2) + 3 + nc1])
    npt <- nc2
    np <- nc1 + nc2
    tfs <- tf(uft)
  }
  
  b <- if (is.null(init)) 
    rep(0, np) else init
  stepf <- 0.5
  
  for (ll in 0:maxiter) {
    z <- .Fortran(
      'crrfsv', as.double(ftime), as.integer(fstatus), as.integer(length(ftime)),
      as.double(cov1), as.integer(np - npt), as.integer(np), as.double(cov2),
      as.integer(npt), as.double(tfs), as.integer(ndf), as.double(uuu),
      as.integer(ncg), as.integer(cengroup), as.double(b), double(1), double(np),
      double(np * np), double(np), double(np), double(np * np),
      PACKAGE = 'cmprsk'
    )
    z <- z[15:17]
    if (max(abs(z[[2]]) * pmax(abs(b), 1)) < max(abs(z[[1]]), 1) * gtol) {
      converge <- TRUE
      break
    }
    
    if (ll == maxiter) {
      converge <- FALSE
      break
    }
    
    h <- z[[3L]]
    dim(h) <- c(np, np)
    sc <- -solve(h, z[[2L]])
    bn <- b + sc
    
    fbn <- .Fortran(
      'crrf', as.double(ftime), as.integer(fstatus), as.integer(length(ftime)),
      as.double(cov1), as.integer(np - npt), as.integer(np), as.double(cov2),
      as.integer(npt), as.double(tfs), as.integer(ndf), as.double(uuu),
      as.integer(ncg), as.integer(cengroup), as.double(bn), double(1), double(np),
      PACKAGE = 'cmprsk'
    )
    fbn <- fbn[[15L]]
    
    i <- 0
    while (is.na(fbn) || fbn > z[[1]] + (1e-04) * sum(sc * z[[2L]])) {
      i <- i + 1
      sc <- sc * stepf
      bn <- b + sc
      fbn <- .Fortran(
        'crrf', as.double(ftime), as.integer(fstatus), as.integer(length(ftime)),
        as.double(cov1), as.integer(np - npt), as.integer(np), as.double(cov2),
        as.integer(npt), as.double(tfs), as.integer(ndf), as.double(uuu),
        as.integer(ncg), as.integer(cengroup), as.double(bn), double(1), double(np),
        PACKAGE = 'cmprsk'
      )
      fbn <- fbn[[15L]]
      
      if (i > 20) 
        break
    }
    
    if (i > 20) {
      converge <- FALSE
      break
    }
    b <- c(bn)
  }
  if (variance) {
    v <- .Fortran(
      'crrvv', as.double(ftime), as.integer(fstatus), as.integer(length(ftime)),
      as.double(cov1), as.integer(np - npt), as.integer(np), as.double(cov2),
      as.integer(npt), as.double(tfs), as.integer(ndf), as.double(uuu),
      as.integer(ncg), as.integer(cengroup), as.double(b), double(np * np),
      double(np * np), double(np * np), double(length(ftime) * (np + 1)),
      double(np), double(np * ncg), double(2 * np), double(ncg * np), integer(ncg),
      double(ncg * np), double(ncg),
      PACKAGE = 'cmprsk'
    )
    v <- v[15:16]
    dim(v[[2L]]) <- dim(v[[1L]]) <- c(np, np)
    h0 <- v[[1L]]
    h <- solve(v[[1L]])
    v <- h %*% v[[2L]] %*% t(h)
    
    r <- .Fortran(
      'crrsr', as.double(ftime), as.integer(fstatus), as.integer(length(ftime)),
      as.double(cov1), as.integer(np - npt), as.integer(np), as.double(cov2),
      as.integer(npt), as.double(tfs), as.integer(ndf), as.double(uuu),
      as.integer(ncg), as.integer(cengroup), as.double(b), double(ndf * np),
      double(np), double(np),
      PACKAGE = 'cmprsk'
    )
    r <- r[[15L]]
    r <- t(matrix(r, nrow = np))
  }
  else {
    v <- h <- h0 <- matrix(NA, np, np)
    r <- NULL
  }
  
  nobs <- length(ftime)
  b0 <- rep(0, length(b))
  
  fb0 <- .Fortran(
    'crrf', as.double(ftime), as.integer(fstatus), as.integer(length(ftime)),
    as.double(cov1), as.integer(np - npt), as.integer(np), as.double(cov2),
    as.integer(npt), as.double(tfs), as.integer(ndf), as.double(uuu), as.integer(ncg),
    as.integer(cengroup), as.double(b0), double(1), double(np),
    PACKAGE = 'cmprsk'
  )
  fb0 <- fb0[[15L]]
  
  bj <- .Fortran(
    'crrfit', as.double(ftime), as.integer(fstatus), as.integer(length(ftime)),
    as.double(cov1), as.integer(np - npt), as.integer(np), as.double(cov2),
    as.integer(npt), as.double(tfs), as.integer(ndf), as.double(uuu), as.integer(ncg),
    as.integer(cengroup), as.double(b), double(ndf), double(np),
    PACKAGE = 'cmprsk'
  )
  bj <- bj[[15L]]
  
  if (nc1 > 0) {
    x1 <- paste(cov1.name, 1:nc1, sep = '')
    cov1.vars <- if (is.null(cov1.vars))
      x1 else ifelse(cov1.vars == '', x1, cov1.vars)
  }
  
  if (nc2 > 0) {
    x1 <- paste(cov2.name, 1:nc2, sep = '')
    cov2.vars <- if (is.null(cov2.vars))
      x1 else ifelse(cov2.vars == '', x1, cov2.vars)
    x1 <- paste('tf', 1:nc2, sep = '')
    x2 <- colnames(tfs)
    if (!is.null(x2)) 
      x1 <- ifelse(x2 == '', x1, x2)
    cov2.vars <- paste(cov2.vars, x1, sep = '*')
  }
  names(b) <- c(cov1.vars, cov2.vars)
  
  structure(
    list(coef = b, loglik = -z[[1L]], score = -z[[2L]], inf = h0,
         var = v, res = r, uftime = uft, bfitj = bj, tfs = as.matrix(tfs),
         converged = converge, call = call, n = nobs, n.missing = nmis,
         loglik.null = -fb0, invinf = h),
    class = 'crr'
  )
}
