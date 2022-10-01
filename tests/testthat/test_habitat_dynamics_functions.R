context('habitat_dynamics_functions-class')

test_that('habitat dynamics functions class works', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k,
                         "fires" = egk_fire)
  
  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                        effect_time = 2)),
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  
})
