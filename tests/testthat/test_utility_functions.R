context('utility-functions')

test_that('utility functions work', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k)
  
  landscape2 <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = "k")
  
  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                 global_stochasticity = 0.05),
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE,
                    demo_stochasticity = "none")
  
  expect_error(simulation(landscape = landscape2,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE,
                    demo_stochasticity = "none"))
  
  
})
