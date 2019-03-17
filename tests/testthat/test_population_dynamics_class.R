context('population_dynamics-class')

test_that('population_dynamics classes work', {

  func <- function(x) x
  
  expect_true(inherits(population_dynamics(func),"population_dynamics"))
  
  expect_true(is.population_dynamics(population_dynamics(func)))

  expect_error(as.population_dynamics(1))
  
  expect_error(population_dynamics(change = NULL,
                                   dispersal = NULL,
                                   modification = NULL,
                                   density_dependence = NULL,
                                   dynamics_order = "change"))
    
  # print(population_dynamics(func))
  
})