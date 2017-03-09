context('terms.inner')

# x <- Surv(time, status(0) == 'event') ~ age + sex
# gsub('Surv\\(.*(?:==|%in%)\\s*[\'\" ]*?(.*)[\'\" ]*\\s*\\)|.', '\\1', deparse(x))
# gsub('Surv\\(.*(?:==|%in%)\\s*[\'\"]?(.*)[\'\"]?\\s*\\)|.', '\\1', deparse(form))

test_that('terms.inner parses formula', {
  
  expect_identical(trimwsq(c(' "x"', 'x\'s', ' x ""', 'x')),
                   c('x', 'x\'s', 'x', 'x'))
  
  l <- list(c('time', 'status', 'x', 'y'), c('0', '1'), c('Censored', 'Event'))
  
  expect_identical(terms.inner(Surv(time, status(0) == 1) ~ x + y), l[-3])
  expect_identical(terms.inner(Surv(time, status(Censored) == Event) ~ x + y), l[-2])
  expect_identical(terms.inner(Surv(time, status("Censored")==Event) ~ x + y), l[-2])
  expect_identical(terms.inner(Surv(time, status(0)%in%'1') ~ x + y), l[-2])
  expect_identical(terms.inner(Surv(time, status()%in%'1') ~ x + y), l[-2])
})
