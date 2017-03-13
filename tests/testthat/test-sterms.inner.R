context('sterms.inner')

test_that('sterms.inner parses Surv formula', {
  
  l <- c('time', 'status', 'x', 'y')
  
  expect_identical(
    sterms.inner(Surv(time, status) ~ x + y),
    l
  )
  
  expect_identical(
    sterms.inner(Surv(time, status == 1) ~ x + y),
    l
  )
  
})
