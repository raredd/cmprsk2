context('cmprsk2 utils')


test_that('binary operators', {
  
  expect_identical(
    NULL %||% 1, 1
  )
  expect_identical(
    1 %||% simpleError(), 1
  )
  
  expect_false(
    islist(data.frame())
  )
  expect_true(
    islist(list())
  )
  
})

test_that('package utils', {
  
  ## trimwsq trims white space and quotes
  expect_identical(
    trimwsq(c(' "x"', 'x\'s', ' x ""', 'x')),
    c('x', 'x\'s', 'x', 'x')
  )
  
  expect_true(
    cmprsk2:::is.loaded('cmprsk2')
  )
  expect_false(
    cmprsk2:::is.loaded('no package name')
  )
  
})
