context('simulation_results-class')

test_that('simulation_results classes work', {
  
  library(raster)
  library(future)
  
  egk_mat_mask <- egk_mat > 0
  
  egk_hab_stack <- stack(replicate(10, egk_hab))
  egk_hab_stack_na <- egk_hab_stack
  egk_hab_stack_na[[5]][1:10] <- NA
  
  disp_bar <- egk_hab
  disp_bar[] <- 0
  disp_bar[cellFromCol(disp_bar,ncol(disp_bar)/2)] <- 1
  
  disp_bar2 <- egk_hab
  disp_bar2[] <- 0
  disp_bar2[sampleRandom(disp_bar2, size = 100, na.rm = TRUE, sp = TRUE)] <- 1
  
  cull <- stack(replicate(3, egk_hab))
  cull[] <- 1
  cull[sample(1:ncell(cull), 50)] <- 0
  
  pop_source <- egk_pop[[3]]
  pop_source[] <- 0
  pop_source[sample(which(getValues(egk_pop[[3]]) >= 2), 3)] <- 1
  #plot(pop_source, box = FALSE, axes = FALSE)
  
  pop_sink <- egk_pop[[3]]
  pop_sink[] <- 0
  pop_sink[sample(which(getValues(egk_pop[[3]]) <= 2),
                  cellStats(pop_source, sum))] <- 1
  #plot(pop_sink, box = FALSE, axes = FALSE)
  
  carrying_cap_fun <- function (landscape, timestep) {
    
    fun <- function(suitability) {
      75 - round(75 * dlogis(suitability, scale = 0.25))
    }
    
    suit <- landscape$suitability
    if (raster::nlayers(suit) > 1) {
      suit <- suit[[timestep]]
    }
    
    calc(suit, fun)
    
  }
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k,
                         fires = egk_fire,
                         source = pop_source,
                         sink = pop_sink,
                         barrier = disp_bar,
                         barrier2 = disp_bar2,
                         cull = cull)
  
  landscape_kfun <- landscape(population = egk_pop,
                              suitability = egk_hab,
                              carrying_capacity = carrying_cap_fun)
  
  landscape_nohab <- landscape(population = egk_pop,
                               suitability = NULL,
                               carrying_capacity = egk_k)
  
  landscape_nohab2 <- landscape(population = egk_pop,
                                suitability = NULL,
                                carrying_capacity = carrying_cap_fun)
  
  landscape_nok <- landscape(population = egk_pop,
                             suitability = egk_hab,
                             carrying_capacity = NULL)
  
  landscape_intk <- landscape(population = egk_pop,
                              suitability = egk_hab,
                              carrying_capacity = 0)
  
  landscape_habstack <- landscape(population = egk_pop,
                                  suitability = egk_hab_stack,
                                  carrying_capacity = egk_k,
                                  fires = egk_fire,
                                  source = pop_source,
                                  sink = pop_sink)
  
  expect_error(landscape(population = egk_pop,
                         suitability = egk_hab_stack_na,
                         carrying_capacity = egk_k,
                         fires = egk_fire,
                         source = pop_source,
                         sink = pop_sink))
  
  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                 dispersal = cellular_automata_dispersal(max_cells = c(0, 20, 0),
                                                                         dispersal_proportion = set_proportion_dispersing(proportions = 0.5),
                                                                         barriers_map = "barrier"),
                                 modification = translocation(source_layer = "source",
                                                              sink_layer = "sink",
                                                              stages = NULL,
                                                              effect_timesteps = 2),
                                 density_dependence = ceiling_density())
  
  pop_dyn2 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = fast_dispersal(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 8000)),
                                  modification = mortality(mortality_layer = 'cull', stages = 2),
                                  density_dependence = NULL)
  
  pop_dyn3 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = kernel_dispersal(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 8000),
                                                               max_distance = 8000,
                                                               arrival_probability = "suitability",
                                                               dispersal_proportion = set_proportion_dispersing(proportions = c(0, 0.5, 1))),
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn3e <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                   dispersal = kernel_dispersal(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 8000),
                                                                max_distance = 8000,
                                                                arrival_probability = "both"),
                                   modification = NULL,
                                   density_dependence = NULL)
  
  pop_dyn4 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = cellular_automata_dispersal(max_cells = c(0, 20, 0),
                                                                          dispersal_proportion = density_dependence_dispersing(),
                                                                          barriers_map = "barrier2"),
                                  modification = translocation(source_layer = "source",
                                                               sink_layer = "sink",
                                                               stages = 3,
                                                               effect_timesteps = 2),
                                  density_dependence = ceiling_density())
  
  pop_dyn4a <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   global_stochasticity = matrix(c(0.00,0.10,0.00,0.00,0.00,0.10,0.10,0.00,0.10), nrow = 3, ncol = 3),
                                                   local_stochasticity = matrix(c(0.00,0.10,0.00,0.00,0.00,0.10,0.10,0.00,0.10), nrow = 3, ncol = 3),
                                                   transition_function = modified_transition(survival_layer = "suitability",
                                                                                             fecundity_layer = "suitability")),
                                   dispersal = cellular_automata_dispersal(max_cells = 20),
                                   modification = NULL,
                                   density_dependence = ceiling_density())
  
  pop_dyn4b <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   transition_function = list(modified_transition(), competition_density())),
                                   dispersal = NULL,
                                   modification = NULL,
                                   density_dependence = NULL)
  
  pop_dyn4c <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   transition_function = competition_density(stages = c(2, 3),
                                                                                             mask = egk_mat_mask)),
                                   dispersal = NULL,
                                   modification = NULL,
                                   density_dependence = NULL)
  
  pop_dyn4d <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                   dispersal = cellular_automata_dispersal(max_cells = 20,
                                                                           carrying_capacity = "suitability"),
                                   modification = NULL,
                                   density_dependence = ceiling_density())
  
  pop_dyn5 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                  transition_function = modified_transition(survival_layer = NULL,
                                                                                            fecundity_layer = NULL)),
                                  dispersal = kernel_dispersal(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 8000),
                                                               max_distance = 8000),
                                  modification = NULL,
                                  density_dependence = ceiling_density(stages = c(2,3)))
  
  pop_dyn6 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = cellular_automata_dispersal(),
                                  modification = translocation(source_layer = "source",
                                                               sink_layer = "sink",
                                                               stages = 3,
                                                               effect_timesteps = 2),
                                  density_dependence = ceiling_density())
  
  pop_dyn7 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = fast_dispersal(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.1,
                                                                                                             normalize = TRUE)),
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn7e <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                   global_stochasticity = matrix(c(0,0,0,0,0,1,1,1,1), nrow = 3, ncol = 3)),
                                   dispersal = fast_dispersal(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 8000)),
                                   modification = NULL,
                                   density_dependence = NULL)
  
  pop_dyn8 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = NULL,
                                  modification = NULL,
                                  density_dependence = NULL)
  
  pop_dyn9 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                  dispersal = NULL,
                                  modification = NULL,
                                  density_dependence = ceiling_density())
  
  pop_dyn10 <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                   dispersal = kernel_dispersal(max_distance = 5000,
                                                                arrival_probability = "carrying_capacity"),
                                   modification = NULL,
                                   density_dependence = NULL)
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn,
                                  habitat_dynamics = NULL,
                                  timesteps = 3)[1],
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn,
                                  habitat_dynamics = NULL,
                                  demo_stochasticity = "deterministic_redistribution",
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(is.simulation_results(simulation(landscape = landscape,
                                               population_dynamics = pop_dyn,
                                               habitat_dynamics = NULL,
                                               demo_stochasticity = "stochastic_redistribution",
                                               timesteps = 3)))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn2,
                                  habitat_dynamics = NULL,
                                  demo_stochasticity = "none",
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn3,
                                  habitat_dynamics = NULL,
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn4,
                                  habitat_dynamics = NULL,
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape_kfun,
                                  population_dynamics = pop_dyn4a,
                                  habitat_dynamics = NULL,
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape_habstack,
                                  population_dynamics = pop_dyn4a,
                                  habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                      effect_time = 2)),
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn4b,
                                  habitat_dynamics = NULL,
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn4c,
                                  habitat_dynamics = NULL,
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn4d,
                                  habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                      effect_time = 2)),
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn5,
                                  habitat_dynamics = NULL,
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn6,
                                  habitat_dynamics = NULL,
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn7,
                                  habitat_dynamics = NULL,
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn,
                                  habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                      effect_time = 2)),
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn2,
                                  habitat_dynamics = list(fire_effects(fire_layers = "fires",
                                                                       effect_time = 2)),
                                  timesteps = 3),
                       "simulation_results"))
  
  expect_true(inherits(simulation(landscape = landscape,
                                  population_dynamics = pop_dyn,
                                  habitat_dynamics = NULL,
                                  timesteps = 3,
                                  replicates = 2),
                       "simulation_results"))
  
  
  expect_error(simulation(landscape = landscape_nohab,
                          population_dynamics = pop_dyn,
                          habitat_dynamics = NULL,
                          timesteps = 3))
  
  expect_error(simulation(landscape = landscape_nok,
                          population_dynamics = pop_dyn,
                          habitat_dynamics = NULL,
                          timesteps = 3))
  
  expect_error(simulation(landscape = landscape_habstack,
                          population_dynamics = pop_dyn2,
                          habitat_dynamics = list(fire_effects(fire_layers = "fires",
                                                               effect_time = 2)),
                          timesteps = 3))
  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn2,
                          habitat_dynamics = list(fire_effects(fire_layers = "fires",
                                                               effect_time = 2)),
                          timesteps = 21))
  
  expect_error(simulation(landscape = landscape_nohab,
                          population_dynamics = pop_dyn3e,
                          habitat_dynamics = NULL,
                          timesteps = 3))
  
  expect_error(simulation(landscape = landscape_nok,
                          population_dynamics = pop_dyn3e,
                          habitat_dynamics = NULL,
                          timesteps = 3))
  
  expect_error(simulation(landscape = landscape_nok,
                          population_dynamics = pop_dyn4,
                          habitat_dynamics = NULL,
                          timesteps = 3))
  
  expect_error(simulation(landscape = landscape_intk,
                          population_dynamics = pop_dyn4,
                          habitat_dynamics = NULL,
                          timesteps = 3))
  
  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn7e,
                          habitat_dynamics = NULL,
                          timesteps = 10))
  
  expect_error(simulation(landscape = landscape_nok,
                          population_dynamics = pop_dyn9,
                          habitat_dynamics = NULL,
                          timesteps = 10))
  
  expect_error(simulation(landscape = landscape_nok,
                          population_dynamics = pop_dyn10,
                          habitat_dynamics = NULL,
                          timesteps = 10))
  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn,
                          habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                              effect_time = 2)),
                          timesteps = 21))
  
  test_simulation <- simulation(landscape = landscape,
                                population_dynamics = pop_dyn,
                                habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                    effect_time = 2)),
                                timesteps = 10,
                                replicates = 3)
  
  plan(multisession)
  test_simulation_par <- simulation(landscape = landscape,
                                population_dynamics = pop_dyn,
                                habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                                    effect_time = 2)),
                                timesteps = 10,
                                replicates = 3)
  
  
  plan(sequential)
  test_simulation2 <- simulation(landscape = landscape,
                                 population_dynamics = pop_dyn8,
                                 habitat_dynamics = NULL,
                                 timesteps = 30,
                                 replicates = 1)
  
  #print(test_simulation)
  
  plot(test_simulation)
  
  plot(test_simulation,
       stage = 0,
       emp = TRUE)
  
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
       type = "raster")
  
  plot(test_simulation[1],
       object = "carrying_capacity",
       type = "raster",
       animate = TRUE)
  
  expect_error(plot(test_simulation,
                    object = "population",
                    type = "raster"))
  
  expect_error(plot(test_simulation[1],
                    object = "population",
                    type = "raster"))
  
  expect_true(inherits(extract_spatial(test_simulation2),
                       "RasterLayer"))
  
  expect_true(inherits(extract_spatial(test_simulation2,
                                       landscape_object = 4),
                       "RasterLayer"))
  
  compare_emp(test_simulation, test_simulation, test_simulation)
  
})
