context('population_dispersal_functions-class')

test_that('population dispersal functions class works', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k,
                         "barriers" = egk_road)
  
  landscape_sing_roads <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k,
                         "barriers" = egk_road[[1]])
  
  landscape_nohab <- landscape(population = egk_pop,
                         suitability = NULL,
                         carrying_capacity = egk_k)
  
  landscape_nok <- landscape(population = egk_pop,
                               suitability = egk_hab,
                               carrying_capacity = NULL)

  pop_dyn_kd <- population_dynamics(change = NULL,
                                 dispersal = kernel_dispersal(exponential_dispersal_kernel(distance_decay = 8000)),
                                 modification = NULL,
                                 density_dependence = NULL)
  
  pop_dyn_kd_cc <- population_dynamics(change = NULL,
                                    dispersal = kernel_dispersal(exponential_dispersal_kernel(distance_decay = 8000),
                                                                 arrival_probability = "carrying_capacity"),
                                    modification = NULL,
                                    density_dependence = NULL)

  pop_dyn_kd_large <- population_dynamics(change = NULL,
                                     dispersal = kernel_dispersal(exponential_dispersal_kernel(distance_decay = 8000), max_distance = 100000),
                                     modification = NULL,
                                     density_dependence = NULL)
  
  pop_dyn_kd_neg <- population_dynamics(change = NULL,
                                     dispersal = kernel_dispersal(exponential_dispersal_kernel(distance_decay = 8000), max_distance = -1),
                                     modification = NULL,
                                     density_dependence = NULL)
  
  pop_dyn_kd_inf <- population_dynamics(change = NULL,
                                     dispersal = kernel_dispersal(exponential_dispersal_kernel(distance_decay = 8000), max_distance = Inf),
                                     modification = NULL,
                                     density_dependence = NULL)
  
  pop_dyn_kd_bad <- population_dynamics(change = NULL,
                                        dispersal = kernel_dispersal(exponential_dispersal_kernel(distance_decay = 8000), max_distance = c(1, 2)),
                                        modification = NULL,
                                        density_dependence = NULL)
 
  pop_dyn_ca <- population_dynamics(change = NULL,
                                    dispersal = cellular_automata_dispersal(),
                                    modification = NULL,
                                    density_dependence = NULL)
  
  pop_dyn_ca_sls <- population_dynamics(change = NULL,
                                        dispersal = cellular_automata_dispersal(min_cells = 0,
                                                                                max_cells = 10),
                                        modification = NULL,
                                        density_dependence = NULL)
  
  pop_dyn_ca_barriers <- population_dynamics(change = NULL,
                                        dispersal = cellular_automata_dispersal(barriers = "barriers",
                                                                                use_suitability = FALSE),
                                        modification = NULL,
                                        density_dependence = NULL)
  
  pop_dyn_ca_minmax <- population_dynamics(change = NULL,
                                    dispersal = cellular_automata_dispersal(min_cells = c(0, 1, 0), max_cells = 2),
                                    modification = NULL,
                                    density_dependence = NULL)
  
  pop_dyn_fd <- population_dynamics(change = NULL,
                                    dispersal = fast_dispersal(),
                                    modification = NULL,
                                    density_dependence = NULL)
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_kd,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_kd_large,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_kd_inf,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn_kd_bad,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
  
  expect_error(simulation(landscape = landscape,
                    population_dynamics = pop_dyn_kd_neg,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE))
  
  expect_error(simulation(landscape = landscape_nohab,
                          population_dynamics = pop_dyn_kd,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
  
  expect_error(simulation(landscape = landscape_nok,
                          population_dynamics = pop_dyn_kd,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
  
  expect_error(simulation(landscape = landscape_nok,
                          population_dynamics = pop_dyn_kd_cc,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
 
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_ca,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)

  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_ca_sls,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_ca_barriers,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim <- simulation(landscape = landscape_sing_roads,
                    population_dynamics = pop_dyn_ca_barriers,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  expect_error(simulation(landscape = landscape_nohab,
                          population_dynamics = pop_dyn_ca,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
  
  expect_error(simulation(landscape = landscape_nok,
                          population_dynamics = pop_dyn_ca,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn_ca_minmax,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
   
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_fd,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)

})
