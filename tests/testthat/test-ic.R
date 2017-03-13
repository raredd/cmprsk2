context('information criterion')

test_that('IC aliases return identical results', {
  
  crr <- crr2(Surv(futime, event(censored) == death) ~ age + sex + abo,
              transplant)[[1L]]
  
  expect_identical(
    AIC(crr), extractAIC(crr)
  )
  expect_identical(
    AIC(crr), extractIC(crr, 'AIC')
  )
  
  expect_identical(
    BIC(crr), extractIC(crr, 'BIC')
  )
  
  expect_identical(
    -2 * crr$loglik, unname(extractIC(crr, p = 0))
  )
  
  expect_identical(
    logLik(crr)[[1L]], crr$loglik
  )
  
})
