context('habitat_dynamics-class')

test_that('habitat_dynamics classes work', {

  func <- function(x) x
  
  expect_true(inherits(as.habitat_dynamics(func),"habitat_dynamics"))

  expect_error(as.habitat_dynamics(1))
    
  print(as.habitat_dynamics(func))

  #expect_error(as.demography(c(1,2,3)))
  
})