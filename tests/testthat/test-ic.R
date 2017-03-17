context('information criterion')

test_that('IC aliases return identical results', {
  
  crr <- crr2(Surv(futime, event(censored) == death) ~ age + sex + abo,
              transplant)[[1L]]
  
  ## ignore names in results
  expect_identical_unnamed <- function(object, expected, ...) {
    testthat::expect_identical(
      unname(object),
      unname(expected),
      ...
    )
  }
  
  expect_identical_unnamed(
    AIC(crr), extractAIC(crr)
  )
  expect_identical_unnamed(
    AIC(crr), extractIC(crr, 'AIC')
  )
  
  expect_identical_unnamed(
    BIC(crr), extractIC(crr, 'BIC')
  )
  
  expect_identical_unnamed(
    -2 * crr$loglik, extractIC(crr, p = 0)
  )
  
  expect_identical_unnamed(
    logLik(crr)[[1L]], crr$loglik
  )
  
})
