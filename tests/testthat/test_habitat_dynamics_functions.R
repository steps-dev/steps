context('habitat_dynamics_functions-class')

test_that('habitat dynamics functions class works', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k,
                         "fires" = egk_fire)
  
  egk_hab_stack <- stack(replicate(10, egk_hab))
  egk_hab_stack_na <- egk_hab_stack
  egk_hab_stack_na[[5]][1:10] <- NA
  
  landscape_stacks <- landscape(population = egk_pop,
                                suitability = egk_hab_stack,
                                carrying_capacity = egk_k,
                                "fires" = egk_fire)
  
  landscape_bad_layers <- landscape(population = egk_pop,
                                suitability = egk_hab_stack,
                                carrying_capacity = egk_k,
                                "fires" = egk_fire[[1]])
  
  landscape_bad_layers2 <- landscape(population = egk_pop,
                                    suitability = egk_hab,
                                    carrying_capacity = egk_k,
                                    "fires" = egk_fire[[1]])
  
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
  
  sim <- simulation(landscape = landscape_stacks,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                        effect_time = 2)),
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  expect_error(simulation(landscape = landscape_bad_layers,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                        effect_time = 2)),
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE))
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = list(fire_effects(fire_layers = "fires",
                                                        effect_time = 2)),
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  expect_error(simulation(landscape = landscape_bad_layers,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = list(fire_effects(fire_layers = "fires",
                                                         effect_time = 2)),
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE))
  
  expect_error(simulation(landscape = landscape_bad_layers2,
                          population_dynamics = pop_dyn,
                          habitat_dynamics = list(fire_effects(fire_layers = "fires",
                                                               effect_time = 2)),
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
  
  
})
