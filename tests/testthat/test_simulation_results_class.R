context('simulation_results-class')

options(mypkg.connection = stdin())

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
  # egk_fire <- stack(dist_list)
  # for (i in c(3,7,14,15,18)) {
  #   egk_fire[[i]][] <- 1
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
                   carrying_capacity = egk_k,
                   fires = egk_fire,
                   source = pop_source,
                   sink = pop_sink)
  
  landscape_nohab <- landscape(population = egk_pop,
                         suitability = NULL,
                         carrying_capacity = egk_k)
  
  landscape_habstack <- landscape(population = egk_pop,
                               suitability = stack(replicate(10, egk_hab)),
                               carrying_capacity = egk_k,
                               fires = egk_fire,
                               source = pop_source,
                               sink = pop_sink)
  
  expect_error(landscape(population = egk_pop,
                         suitability = NULL,
                         carrying_capacity = function(k) k*10))

  landscape_kfun <- landscape(population = egk_pop,
                               suitability = egk_hab,
                               carrying_capacity = function(k) k*10)

  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                 dispersal = cellular_automata_dispersal(dispersal_distance=c(0, 10, 0),
                                                                    dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                                    dispersal_proportion=c(0, 0.25, 0),
                                                                    barrier_type = 1,
                                                                    use_barriers = TRUE,
                                                                    barriers_map = disp.bar),
                                 modification = translocation(source_layer = "source",
                                                     sink_layer = "sink",
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
                                  dispersal = cellular_automata_dispersal(dispersal_distance=c(0, 10, 0),
                                                                     dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                                     dispersal_proportion=c(0, 0.25, 0),
                                                                     barrier_type = 0,
                                                                     use_barriers = TRUE,
                                                                     barriers_map = disp.bar2),
                                  modification = translocation(source_layer = "source",
                                                      sink_layer = "sink",
                                                      stages = 3,
                                                      effect_timesteps = 2),
                                  density_dependence = population_cap())
  
  pop_dyn4a <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   global_stochasticity = matrix(c(0.00,0.10,0.00,0.00,0.00,0.10,0.10,0.00,0.10), nrow = 3, ncol = 3),
                                                   local_stochasticity = matrix(c(0.00,0.10,0.00,0.00,0.00,0.10,0.10,0.00,0.10), nrow = 3, ncol = 3),
                                                   transition_function = modified_transition(egk_mat,
                                                                       survival_layer = "suitability",
                                                                       fecundity_layer = "suitability")),
                                  dispersal = cellular_automata_dispersal(dispersal_distance=c(10),
                                                                          dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                                          dispersal_proportion=c(0.25)),
                                  modification = NULL,
                                  density_dependence = population_cap())
  
  pop_dyn5 <- population_dynamics(change = growth(transition_matrix = egk_mat, demographic_stochasticity = FALSE),
                                  dispersal = kernel_dispersal(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.1)),
                                  modification = NULL,
                                  density_dependence = population_cap(stages = c(2,3)))

  pop_dyn6 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = cellular_automata_dispersal(dispersal_distance=list(0, 10, 0),
                                                                          dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                                          dispersal_proportion=list(0, 0.25, 0)),
                                  modification = translocation(source_layer = "source",
                                                               sink_layer = "sink",
                                                               stages = 3,
                                                               effect_timesteps = 2),
                                  density_dependence = population_cap())
  
  pop_dyn7 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1,
                                                                                                           normalize = TRUE)),
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn7e <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   global_stochasticity = matrix(c(0,0,0,0), nrow = 2, ncol = 2)),
                                  dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1)),
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn7e2 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   local_stochasticity = matrix(c(0,0,0,0), nrow = 2, ncol = 2)),
                                   dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1)),
                                   modification = NULL,
                                   density_dependence = NULL)
  
  pop_dyn7e3 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   global_stochasticity = matrix(c(0,0,0,0,0,1,1,1,1), nrow = 3, ncol = 3)),
                                   dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1)),
                                   modification = NULL,
                                   density_dependence = NULL)
  
  pop_dyn7e4 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                    local_stochasticity = matrix(c(0,0,0,0,0,1,1,1,1), nrow = 3, ncol = 3)),
                                    dispersal = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1)),
                                    modification = NULL,
                                    density_dependence = NULL)
  
  pop_dyn8 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = NULL,
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
  
  expect_true(inherits(simulation(landscape = landscape_kfun,
                                  population_dynamics = pop_dyn4a,
                                  habitat_dynamics = NULL,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape_habstack,
                                  population_dynamics = pop_dyn4a,
                                  habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                      effect_time = 2)),
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
                                  habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                      effect_time = 2)),
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn2,
                                  habitat_dynamics = list(fire_effects(fire_layers = "fires",
                                                                      lag = 2)),
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
                                  habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                      effect_time = 2)),
                                  timesteps = 30),
                       "simulation_results"))
   
  test_simulation <- simulation(landscape = landscape,
                                population_dynamics = pop_dyn,
                                habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                    effect_time = 2)),
                                timesteps = 10,
                                replicates = 3)
  
  test_simulation2 <- simulation(landscape = landscape,
                                population_dynamics = pop_dyn8,
                                habitat_dynamics = NULL,
                                timesteps = 30,
                                replicates = 1)
    
  print(test_simulation)

  plot(test_simulation)
  
  plot(test_simulation,
       stage = 0)
  
  plot(test_simulation,
       stage = 2)
  
  plot(test_simulation[1])
  
  plot(test_simulation[c(2:3)])
  
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
  
  # f <- file()
  # options(mypkg.connection = f)
  # ans <- paste(c("yes", "no"), collapse = "\n")
  # write(ans, f)
  # 
  # plot(test_simulation2,
  #                   object = "population",
  #                   type = "raster",
  #                   stage = 0,
  #                   timesteps = 1:30)
  # 
  # options(mypkg.connection = stdin())
  # close(f)

  expect_true(inherits(extract_spatial(test_simulation2),
                       "RasterLayer"))
  
  expect_true(inherits(extract_spatial(test_simulation2,
                                       landscape_object = 4),
                       "RasterLayer"))
  
  print(exponential_dispersal_kernel())
  
  print(modified_transition(egk_mat))
  
})
