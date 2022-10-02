context('simulation_results-class')

test_that('simulation_results classes work', {
  
  library(raster)
  library(future)

  egk_hab_stack <- stack(replicate(10, egk_hab))
  egk_hab_stack_na <- egk_hab_stack
  egk_hab_stack_na[[5]][1:10] <- NA

  egk_pop <- round(egk_pop * .75)

  carrying_cap_fun <- function (landscape, timestep) {
    
    fun <- function(suitability) {
      75 - round(75 * dlogis(suitability, scale = 0.25))
    }
    
    suit <- landscape$suitability
    if (raster::nlayers(suit) > 1) {
      suit <- suit[[timestep]]
    }
    
    raster::calc(suit, fun)
    
  }

  # Rmax_pop <- abs(eigen(egk_mat)$values[1])
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k,
                         "fires" = egk_fire)
  
  landscape_kfun <- landscape(population = egk_pop,
                              suitability = egk_hab,
                              carrying_capacity = carrying_cap_fun,
                              "fires" = egk_fire)

  landscape_badk <- landscape(population = egk_pop,
                              suitability = egk_hab,
                              carrying_capacity = "0")
  
  landscape_stacks <- landscape(population = egk_pop,
                                suitability = egk_hab_stack,
                                carrying_capacity = egk_k,
                                fires = egk_fire)

  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL)

  pop_dyn2 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   transition_function = competition_density(R_max = Rmax_pop)),
                                   dispersal = NULL,
                                   modification = NULL,
                                   density_dependence = NULL)

  expect_error(landscape(population = egk_pop,
                         suitability = egk_hab_stack_na,
                         carrying_capacity = egk_k,
                         fires = egk_fire))

  expect_error(simulation(landscape = landscape_badk,
                          population_dynamics = pop_dyn,
                          habitat_dynamics = NULL,
                          timesteps = 3,
                          verbose = FALSE))
  
  expect_error(simulation(landscape = landscape_stacks,
                          population_dynamics = pop_dyn,
                          habitat_dynamics = list(fire_effects(fire_layers = "fires",
                                                               effect_time = 2)),
                          timesteps = 3,
                          verbose = FALSE))
  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn,
                          habitat_dynamics = list(fire_effects(fire_layers = "fires",
                                                               effect_time = 2)),
                          timesteps = 21,
                          verbose = FALSE))

  test_simulation <- simulation(landscape = landscape_kfun,
                                population_dynamics = pop_dyn,
                                habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                    effect_time = 2)),
                                timesteps = 10,
                                replicates = 3,
                                verbose = TRUE)
  
  plan(multisession)
  test_simulation_par <- simulation(landscape = landscape,
                                population_dynamics = pop_dyn,
                                habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                    effect_time = 2)),
                                timesteps = 10,
                                replicates = 2,
                                verbose = FALSE)

  plan(multisession)  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn,
                          habitat_dynamics = NULL,
                          timesteps = 21,
                          verbose = FALSE))
  
  plan(multisession)  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn2,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 2,
                          verbose = FALSE))

  plan(sequential)
  
  plot(test_simulation, replicates = 1:3)

  expect_true(inherits(extract_spatial(test_simulation),
                       "RasterLayer"))
  
  expect_true(inherits(extract_spatial(test_simulation,
                                       landscape_object = "fires"),
                       "RasterLayer"))
  
  expect_true(inherits(extract_spatial(test_simulation,
                                       landscape_object = "suitability"),
                       "RasterLayer"))
  
  expect_true(inherits(extract_spatial(test_simulation,
                                       landscape_object = 4),
                       "RasterLayer"))
  
  steps:::is.simulation_results(test_simulation)

})
