context('population_dynamics-class')

test_that('population_dynamics classes work', {

  func <- function(x) x
  
  expect_true(inherits(as.population_dynamics(func),"population_dynamics"))
  
  expect_true(is.population_dynamics(as.population_dynamics(func)))

  expect_error(as.population_dynamics(1))
    
  print(as.population_dynamics(func))

  #expect_error(as.population(c(1,2,3)))
  
})