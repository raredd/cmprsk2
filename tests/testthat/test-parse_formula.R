context('crr2/cuminc2 formula parsing')


test_that('crr2/cuminc2 formulas are properly parsed', {
  
  data('transplant', package = 'survival')
  tp <- within(transplant, {
    futime_bad <- futime
    futime_bad[1L] <- -futime_bad[1L]
  })
  
  f0 <- Surv(futime, event) ~ 1
  f1 <- Surv(futime, event == 0) ~ 1
  f2 <- Surv(futime, event(0)) ~ sex
  f3 <- Surv(futime, event(censored)) ~ sex
  f4 <- Surv(futime, event(censor) == deaths) ~ sex + abo
  f5 <- Surv(futime, event(censored) == death) ~ sex + strata(abo)
  f6 <- Surv(futime_bad, event(censored)) ~ sex
  
  
  ## check for cencode in formula
  expect_error(
    parse_formula(f0),
    regexp = 'Invalid formula'
  )
  
  expect_error(
    parse_formula(f1),
    regexp = 'Invalid formula'
  )
  
  
  ## cencode and failcode
  expect_identical(
    parse_formula(f2)$cencode,
    '0'
  )
  
  expect_identical(
    parse_formula(f2)$failcode,
    NULL
  )
  
  expect_identical(
    c(parse_formula(f4)$cencode,
      parse_formula(f4)$failcode),
    c('censor', 'deaths')
  )
  
  
  ## rhs and strata
  expect_identical(
    parse_formula(f2)$rhs,
    'sex'
  )
  
  expect_identical(
    parse_formula(f4)$rhs,
    c('sex', 'abo')
  )
  
  expect_identical(
    parse_formula(f5)$rhs,
    'sex'
  )
  
  expect_identical(
    parse_formula(f5)$strata,
    'abo'
  )
  
  
  ## data
  expect_warning(
    parse_formula(f2, tp),
    regexp = '0.*not found.*tp.*event'
  )
  
  expect_warning(
    parse_formula(f4, tp),
    regexp = 'censor.*deaths.*not found.*tp.*event'
  )
  
  expect_error(
    parse_formula(f6, tp),
    regexp = 'numeric values >= 0'
  )
  
})
