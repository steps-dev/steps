context('demography_dynamics-class')

test_that('demography_dynamics classes work', {

  func <- function(x) x
  
  expect_true(inherits(as.demography_dynamics(func),"demography_dynamics"))

  expect_error(as.demography_dynamics(1))
    
  print(as.demography_dynamics(func))

  #expect_error(as.demography(c(1,2,3)))
  
})