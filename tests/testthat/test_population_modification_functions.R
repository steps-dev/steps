context('population_modification_functions-class')

test_that('population modification functions class works', {
  
  library(raster)
  library(future)
  
  pop_origin <- egk_pop[[3]]
  pop_origin[] <- 0
  pop_origin[sample(which(getValues(egk_pop[[3]]) >= 2), 3)] <- 1
  
  pop_destination <- egk_pop[[3]]
  pop_destination[] <- 0
  pop_destination[sample(which(getValues(egk_pop[[3]]) <= 2),
                         cellStats(pop_origin, sum))] <- 1
  
  cull <- stack(replicate(20, egk_hab))
  cull[] <- 1
  cull[sample(1:ncell(cull), 50)] <- 0
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k,
                         "origin" = pop_origin,
                         "destination" = pop_destination,
                         "cull" = cull)
  
  pop_dyn_trans <- population_dynamics(change = NULL,
                                 dispersal = NULL,
                                 modification = translocation(origins_layer = "origin",
                                                              destinations_layer = "destination",
                                                              stages = NULL,
                                                              effect_timesteps = c(2, 5, 8)),
                                 density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_trans,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  pop_dyn_mort <- population_dynamics(change = NULL,
                                    dispersal = NULL,
                                    modification = mortality(mortality_layer = 'cull', stages = 2),
                                    density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_mort,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)

})
