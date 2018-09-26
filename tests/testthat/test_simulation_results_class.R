context('simulation_results-class')

test_that('simulation_results classes work', {
  
  library(raster)
  library(rgdal)
  library(fields)

  disp.bar <- egk_hab
  disp.bar[] <- 0
  disp.bar[cellFromCol(disp.bar,ncol(disp.bar)/2)] <- 1
  
  disp.bar2 <- egk_hab
  disp.bar2[] <- 0
  disp.bar2[sampleRandom(disp.bar2, size = 100, na.rm = TRUE, sp = TRUE)] <- 1
  
  # dist_list <- list()
  # for (i in 1:20) {
  #   r2 <- egk_hab
  #   r2[] <- 1
  #   cells <- sample(c(1:ncell(r2)), sample(5:10, 1))
  #   fireprob <- abs(rnorm(16, 0, 1))
  #   fireprob <- fireprob/max(fireprob)
  #   r2[c(adjacent(egk_hab, cells, directions = 16, pairs = FALSE), cells)]  <- fireprob
  #   dist_list[[i]] <- r2
  # }
  # egk_dist <- stack(dist_list)
  # for (i in c(3,7,14,15,18)) {
  #   egk_dist[[i]][] <- 1
  # }


  # surv_fec <- list(dist.s, dist.s, dist.s)
  # surv_fec2 <- list(dist.s, dist.s, NULL)
  
  pop_source <- egk_pop[[3]]
  pop_source[] <- 0
  pop_source[sample(which(getValues(egk_pop[[3]]) >= 2), 3)] <- 1
  #plot(pop_source, box = FALSE, axes = FALSE)
  
  pop_sink <- egk_pop[[3]]
  pop_sink[] <- 0
  pop_sink[sample(which(getValues(egk_pop[[3]]) <= 2),
                        cellStats(pop_source, sum))] <- 1
  #plot(pop_sink, box = FALSE, axes = FALSE)

  landscape <- landscape(population = egk_pop,
                   suitability = egk_hab,
                   carrying_capacity = egk_k)
  
  landscape_nohab <- landscape(population = egk_pop,
                         suitability = NULL,
                         carrying_capacity = egk_k)

  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                 dispersal = cellular_automata_dispersal(dispersal_distance=list(0, 10, 0),
                                                                    dispersal_kernel=list(0, exp(-c(0:9)^1/3.36), 0),
                                                                    dispersal_proportion=list(0, 0.35*0.714, 0),
                                                                    barrier_type = 1,
                                                                    use_barriers = TRUE,
                                                                    barriers_map = disp.bar),
                                 modification = translocation(source_layer = pop_source,
                                                     sink_layer = pop_sink,
                                                     stages = NULL,
                                                     effect_timesteps = 2),
                                 density_dependence = population_cap())
  
  pop_dyn2 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1)),
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn3 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = kernel_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                          arrival_probability="suitability",
                                                          demographic_stochasticity = TRUE),
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn3e <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = kernel_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                               arrival_probability="both"),
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn4 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = cellular_automata_dispersal(dispersal_distance=list(0, 10, 0),
                                                                     dispersal_kernel=list(0, exp(-c(0:9)^1/3.36), 0),
                                                                     dispersal_proportion=list(0, 0.35*0.714, 0),
                                                                     barrier_type = 1,
                                                                     use_barriers = TRUE,
                                                                     barriers_map = disp.bar2),
                                  modification = translocation(source_layer = pop_source,
                                                      sink_layer = pop_sink,
                                                      stages = 3,
                                                      effect_timesteps = 2),
                                  density_dependence = population_cap())
  
  pop_dyn5 <- population_dynamics(change = growth(transition_matrix = egk_mat, demographic_stochasticity = FALSE),
                                  dispersal = kernel_dispersal(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.1),
                                                          stages = c(2,3)),
                                  modification = NULL,
                                  density_dependence = population_cap(stages = c(2,3)))

  pop_dyn6 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = cellular_automata_dispersal(dispersal_distance=list(0, 10, 0),
                                                                          dispersal_kernel=list(0, exp(-c(0:9)^1/3.36), 0),
                                                                          dispersal_proportion=list(0, 0.35*0.714, 0)),
                                  modification = translocation(source_layer = pop_source,
                                                               sink_layer = pop_sink,
                                                               stages = 3,
                                                               effect_timesteps = 2),
                                  density_dependence = population_cap())
  
  pop_dyn7 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1,
                                                                                                           normalize = TRUE),
                                                             stages = c(2,3)),
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn7e <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   global_stochasticity = matrix(c(0,0,0,0), nrow = 2, ncol = 2)),
                                  dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                             stages = c(2,3)),
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn7e2 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   local_stochasticity = matrix(c(0,0,0,0), nrow = 2, ncol = 2)),
                                   dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                              stages = c(2,3)),
                                   modification = NULL,
                                   density_dependence = NULL)
  
  pop_dyn7e3 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   global_stochasticity = matrix(c(0,0,0,0,0,1,1,1,1), nrow = 3, ncol = 3)),
                                   dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                              stages = c(2,3)),
                                   modification = NULL,
                                   density_dependence = NULL)
  
  pop_dyn7e4 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                    local_stochasticity = matrix(c(0,0,0,0,0,1,1,1,1), nrow = 3, ncol = 3)),
                                    dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                               stages = c(2,3)),
                                    modification = NULL,
                                    density_dependence = NULL)
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn,
                                  habitat_dynamics = NULL,
                                  timesteps = 5)[1],
                       "simulation_results"))
    
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn,
                                  habitat_dynamics = NULL,
                                  timesteps = 5),
                       "simulation_results"))
  
  expect_true(is.simulation_results(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn,
                                  habitat_dynamics = NULL,
                                  timesteps = 10)))

  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn2,
                                  habitat_dynamics = NULL,
                                  timesteps = 10),
                       "simulation_results"))
    
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn3,
                                  habitat_dynamics = NULL,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn4,
                                  habitat_dynamics = NULL,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn5,
                                  habitat_dynamics = NULL,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn6,
                                  habitat_dynamics = NULL,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn7,
                                  habitat_dynamics = NULL,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn,
                                  habitat_dynamics = list(disturbance_fires(habitat_suitability = egk_hab,
                                                                       disturbance_layers = egk_dist,
                                                                       effect_time = 2)),
                                  timesteps = 10),
                       "simulation_results"))

  expect_true(inherits(simulation(landscape = landscape,
                                population_dynamics = pop_dyn,
                                habitat_dynamics = NULL,
                                timesteps = 10,
                                replicates = 2),
                       "simulation_results"))
  
  expect_error(inherits(simulation(landscape = landscape_nohab,
                                  population_dynamics = pop_dyn3e,
                                  habitat_dynamics = NULL,
                                  timesteps = 10),
                       "simulation_results"))

  expect_error(inherits(simulation(landscape = landscape,
                                   population_dynamics = pop_dyn7e,
                                   habitat_dynamics = NULL,
                                   timesteps = 10),
                        "simulation_results"))
  
  expect_error(inherits(simulation(landscape = landscape,
                                   population_dynamics = pop_dyn7e2,
                                   habitat_dynamics = NULL,
                                   timesteps = 10),
                        "simulation_results"))
  
  expect_error(inherits(simulation(landscape = landscape,
                                   population_dynamics = pop_dyn7e3,
                                   habitat_dynamics = NULL,
                                   timesteps = 10),
                        "simulation_results"))
  
  expect_error(inherits(simulation(landscape = landscape,
                                   population_dynamics = pop_dyn7e4,
                                   habitat_dynamics = NULL,
                                   timesteps = 10),
                        "simulation_results"))
  
  expect_error(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn,
                                  habitat_dynamics = list(disturbance_fires(habitat_suitability = egk_hab,
                                                                            disturbance_layers = egk_dist,
                                                                            effect_time = 2)),
                                  timesteps = 30),
                       "simulation_results"))
   
  test_simulation <- simulation(landscape = landscape,
                                population_dynamics = pop_dyn,
                                habitat_dynamics = list(disturbance_fires(habitat_suitability = egk_hab,
                                                                          disturbance_layers = egk_dist,
                                                                          effect_time = 2)),
                                timesteps = 10,
                                replicates = 5)
    
  print(test_simulation)

  plot(test_simulation)
  
  plot(test_simulation,
       stage = 0)
  
  plot(test_simulation,
       stage = 2)
  
  plot(test_simulation[1])
  
  plot(test_simulation[c(2:5)])
  
  plot(test_simulation[1],
       object = "population",
       type = "raster",
       stage = 2)
  
  plot(test_simulation[1],
       type = "raster",
       stage = 2,
       animate = TRUE)
  
  plot(test_simulation[1],
       type = "raster",
       stage = 0,
       animate = TRUE)
  
  plot(test_simulation[1],
       object = "population",
       type = "graph",
       stage = 0)
  
  plot(test_simulation[1],
       object = "population",
       type = "raster",
       stage = 0)
  
  plot(test_simulation[1],
       object = "population",
       type = "graph",
       stage = 2)
  
  plot(test_simulation[1],
       object = "suitability")
  
  plot(test_simulation[1],
       object = "suitability",
       animate = TRUE)
  
  plot(test_simulation[1],
       object = "carrying_capacity")
  
  plot(test_simulation[1],
       object = "carrying_capacity",
       animate = TRUE)
    
  expect_error(plot(test_simulation,
                    object = "population",
                    type = "raster"))
  
  expect_error(plot(test_simulation,
                    object = "suitability"))
  
  expect_error(plot(test_simulation,
                    object = "carrying_capacity"))
  
  expect_error(plot(test_simulation[1],
                    object = "population",
                    type = "raster"))

})
