context('demography-class')

test_that('demography classes work', {
  
  # the types of demography
  mat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
  tmat <- as.demography(mat)
  
  # check as.demography won't handle a silly function
  expect_error(as.demography(mat[,-1]))
  expect_error(as.demography(mat[-1,]))
  
  # check they have the right class
  expect_s3_class(tmat, 'demography')
  
  # check is.demography works on demographys
  expect_true(is.demography(tmat))
  
})