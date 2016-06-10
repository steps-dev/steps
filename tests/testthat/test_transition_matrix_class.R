context('transition_matrix-class')

test_that('transition_matrix classes work', {
  
  # the types of transition_matrix
  mat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
  tmat <- as.transition_matrix(mat)
  
  # check as.transition_matrix won't handle a silly function
  expect_error(as.transition_matrix(mat[,-1]))
  expect_error(as.transition_matrix(mat[-1,]))
  
  # check they have the right class
  expect_s3_class(tmat, 'transition_matrix')
  
  # check is.transition_matrix works on transition_matrixs
  expect_true(is.transition_matrix(tmat))
  
  # check is.transition_matrix works on non-transition_matrixs
  expect_error(is.transition_matrix(list()))
  expect_error(is.transition_matrix(NA))
  expect_error(is.transition_matrix(NULL))
})