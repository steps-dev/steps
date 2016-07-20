context('transition-class')

test_that('transition classes work', {
  
  # the types of transition
  mat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
  tmat <- as.transition(mat)
  
  # check as.transition won't handle a silly function
  expect_error(as.transition(mat[,-1]))
  expect_error(as.transition(mat[-1,]))
  
  # check they have the right class
  expect_s3_class(tmat, 'transition')
  
  # check is.transition works on transitions
  expect_true(is.transition(tmat))
  
  # check is.transition works on non-transitions
  expect_error(is.transition(list()))
  expect_error(is.transition(NA))
  expect_error(is.transition(NULL))
})