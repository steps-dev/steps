context('dispersal_kernel_functions-class')

test_that('dispersal kernel functions class works', {
  
  library(terra)
  library(future)
  
  landscape <- landscape(population = egk_pop(),
                         suitability = egk_hab(),
                         carrying_capacity = egk_k())
  
  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                 dispersal = kernel_dispersal(exponential_dispersal_kernel(distance_decay = 8000)),
                                 modification = NULL,
                                 density_dependence = NULL)
  
  pop_dyn_norm <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                 dispersal = kernel_dispersal(exponential_dispersal_kernel(distance_decay = 8000,
                                                                                           normalize = TRUE)),
                                 modification = NULL,
                                 density_dependence = NULL)
  
  expect_no_error(sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE))
  
  expect_no_error(sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_norm,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE))
  
  
})
