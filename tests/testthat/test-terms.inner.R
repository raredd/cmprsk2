context('terms.inner')

test_that('terms.inner parses crr2 formula', {
  
  l <- list(c('time', 'status', 'x', 'y'), c('0', '1'), c('Censored', 'Event'))
  attr(l, 'dots') <- FALSE
  
  ## terms.inner adds attributes not needed in tests
  expect_identical_no_attr <- function(object, expected, ...) {
    attributes(object) <- attributes(expected) <- NULL
    expect_identical(object, expected, ...)
  }
  
  expect_identical_no_attr(
    terms.inner(Surv(time, status(0) == 1) ~ x + y),
    l[-3L]
  )
  
  expect_identical_no_attr(
    terms.inner(Surv(time, status(Censored) == Event) ~ x + y),
    l[-2L]
  )
  
  expect_identical_no_attr(
    terms.inner(Surv(time, status("Censored")==Event) ~ x + y),
    l[-2L]
  )
  
  expect_identical_no_attr(
    terms.inner(Surv(time, status(0)%in%'1') ~ x + y),
    l[-3L]
  )
  
  expect_identical_no_attr(
    terms.inner(Surv(time, status('0')%in%'1') ~ x + y),
    l[-3L]
  )
  
})
