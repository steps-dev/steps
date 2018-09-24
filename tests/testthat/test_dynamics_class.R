context('dynamics-class')

test_that('dynamics classes work', {
  library(raster)
  
  func <- function(x) x
  
  expect_true(inherits(dynamics(as.habitat_dynamics(func),
                                      as.population_dynamics(func),
                                      as.demography_dynamics(func)
                                      ),
                       "dynamics")
              )
  
  expect_true(is.dynamics(dynamics(as.habitat_dynamics(func),
                                         as.population_dynamics(func),
                                         as.demography_dynamics(func)
                                         )
                          )
              )
  
  expect_error(inherits(dynamics(as.habitat_dynamics(func),
                                      as.population_dynamics(func)
                                      ),
                       "dynamics")
              )
 
  expect_error(dynamics(as.habitat_dynamics(func),
                                      as.population_dynamics(func),
                                      as.demography_dynamics(func),
                                      order = "one"
                                      )
               ) 
  
  
  print(dynamics(as.habitat_dynamics(func),
                      as.population_dynamics(func),
                      as.demography_dynamics(func)
                      )
        )

})
 