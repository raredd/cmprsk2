context('cuminc2 formula')


test_that('cuminc2 formula is properly recognized', {
  
  ## crr2 formula
  
  expect_false(
    is.cuminc2(
      Surv(t, s) ~ 1
    )
  )
  
  expect_false(
    is.cuminc2(
      Surv(t, s) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0)) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, I(s %in% '1')) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0)==1) ~ 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0) == "1") ~ 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0) == 1) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0)%in%1) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0)%in%"death") ~ factor(1) + 1
    )
  )
  
  
  ## cuminc2 formula
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0)) ~ 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0)) ~ 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0)) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0)) ~ factor(1) + 1
    )
  )
  
  expect_true(
    is.cuminc2(
      Surv(t, s(0) ) ~ factor(1) + 1
    )
  )
})
