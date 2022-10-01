context('population_dispersal_functions-class')

test_that('population dispersal functions class works', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k)
  
  pop_dyn_kd <- population_dynamics(change = NULL,
                                 dispersal = kernel_dispersal(exponential_dispersal_kernel(distance_decay = 8000)),
                                 modification = NULL,
                                 density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_kd,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  pop_dyn_ca <- population_dynamics(change = NULL,
                                    dispersal = cellular_automata_dispersal(),
                                    modification = NULL,
                                    density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_ca,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  pop_dyn_fd <- population_dynamics(change = NULL,
                                    dispersal = fast_dispersal(),
                                    modification = NULL,
                                    density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_fd,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
})
