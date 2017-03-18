context('crr2 formula')


test_that('crr2 formula is properly recognized', {
  
  expect_false(
    is.crr2(
      Surv(t, s) ~ 1
    )
  )
  
  expect_false(
    is.crr2(
      Surv(t, s) ~ factor(1) + 1
    )
  )
  
  expect_false(
    is.crr2(
      Surv(t, s(0)) ~ factor(1) + 1
    )
  )
  
  expect_false(
    is.crr2(
      Surv(t, I(s %in% '1')) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.crr2(
      Surv(t, s(0)==1) ~ 1
    )
  )
  
  expect_true(
    is.crr2(
      Surv(t, s(0) == "1") ~ 1
    )
  )
  
  expect_true(
    is.crr2(
      Surv(t, s(0) == 1) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.crr2(
      Surv(t, s(0)%in%1) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.crr2(
      Surv(t, s(0)%in%"death") ~ factor(1) + 1
    )
  )
  
})
