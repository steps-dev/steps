context('simulation_results-class')

test_that('simulation_results classes work', {
  
  library(raster)
  library(rgdal)
  library(fields)

  disp.bar <- egk_hab*0
  disp.bar[cellFromCol(disp.bar,ncol(disp.bar)/2)] <- 1
  disp.bar2 <- egk_hab*0
  disp.bar2[sampleRandom(disp.bar2, size=100, na.rm=TRUE, sp=TRUE)] <- 1
  
  dist_list <- list()
  
  for (i in 1:10) {
    
    r2 <- egk_hab
    r2[] <- 1
    cells <- sample(c(1:ncell(r2)), 10)
    r2[c(adjacent(egk_hab, cells, directions=8, pairs=FALSE),cells)]  <- 0.5
    dist_list[[i]] <- r2
    
  }
  
  dist.s <- stack(dist_list)

  surv_fec <- list(dist.s, dist.s, dist.s)
  surv_fec2 <- list(dist.s, dist.s, NULL)
  
  pop_source <- egk_pop[[3]]
  pop_source[] <- 0
  pop_source[sample(which(getValues(egk_pop[[3]]) >= 2), 3)] <- 1
  #plot(pop_source, box = FALSE, axes = FALSE)
  
  pop_sink <- egk_pop[[3]]
  pop_sink[] <- 0
  pop_sink[sample(which(getValues(egk_pop[[3]]) <= 2),
                        cellStats(pop_source, sum))] <- 1
  #plot(pop_sink, box = FALSE, axes = FALSE)

  b_hab <- habitat(habitat_suitability = egk_hab,
                         carrying_capacity = egk_k)
  
  b_hab2 <- habitat(habitat_suitability = egk_hab,
                          carrying_capacity = NULL)
  
  b_pop <- population(egk_pop)
 
  b_dem <- demography(transition_matrix = egk_mat,
                            scale = "local",
                            habitat_suitability = egk_hab)
 
  b_dem2 <- demography(transition_matrix = egk_mat)
  
  b_state <- state(population = b_pop,
                   habitat = b_hab,
                   demography = b_dem)
  
  b_state2 <- state(population = b_pop,
                    habitat = b_hab,
                    demography = b_dem2)
  
  b_state3 <- state(population = b_pop,
                    habitat = b_hab2,
                    demography = b_dem)
  
  hab_dyn <- habitat_dynamics(disturbance_fires(habitat_suitability = egk_hab,
                                                  disturbance_layers = dist.s,
                                                  effect_time = 2))
  
  dem_dyn <- demography_dynamics(environmental_stochasticity(transition_matrix = egk_mat,
                                                             stochasticity = 0.5),
                                 density_dependence(transition_matrix = egk_mat,
                                                    fecundity_fraction = 0.8,
                                                    survival_fraction = 0.8))
    
  dem_dyn2 <- demography_dynamics(surv_fec_modify(transition_matrix = egk_mat,
                                                  surv_layers = surv_fec,
                                                  fec_layers = surv_fec))
  
  dem_dyn3 <- demography_dynamics(surv_fec_modify(transition_matrix = egk_mat,
                                                  surv_layers = surv_fec2,
                                                  fec_layers = surv_fec))

  pop_dyn <- population_dynamics(change = simple_growth(),
                                 disp = cellular_automata_dispersal(dispersal_distance=list(0, 10, 0),
                                                                    dispersal_kernel=list(0, exp(-c(0:9)^1/3.36), 0),
                                                                    dispersal_proportion=list(0, 0.35*0.714, 0),
                                                                    barrier_type = 1,
                                                                    use_barriers = TRUE,
                                                                    barriers_map = disp.bar),
                                 mod = translocation(source_layer = pop_source,
                                                     sink_layer = pop_sink,
                                                     stages = NULL,
                                                     effect_timesteps = 2),
                                 dens_dep = ceiling_density_dependence())
  
  pop_dyn2 <- population_dynamics(change = simple_growth(demo_stoch = TRUE),
                                  disp = fast_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1)),
                                  mod = NULL,
                                  dens_dep = NULL)
  
  pop_dyn3 <- population_dynamics(change = simple_growth(demo_stoch = TRUE),
                                  disp = kernel_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                          arrival_probability="habitat_suitability"),
                                  mod = NULL,
                                  dens_dep = NULL)
  
  pop_dyn4 <- population_dynamics(change = simple_growth(),
                                  disp = cellular_automata_dispersal(dispersal_distance=list(0, 10, 0),
                                                                     dispersal_kernel=list(0, exp(-c(0:9)^1/3.36), 0),
                                                                     dispersal_proportion=list(0, 0.35*0.714, 0),
                                                                     barrier_type = 1,
                                                                     use_barriers = TRUE,
                                                                     barriers_map = disp.bar2),
                                  mod = translocation(source_layer = pop_source,
                                                      sink_layer = pop_sink,
                                                      stages = 3,
                                                      effect_timesteps = 2),
                                  dens_dep = ceiling_density_dependence())
  
  pop_dyn5 <- population_dynamics(change = simple_growth(demo_stoch = TRUE),
                                  disp = kernel_dispersal(dispersal_kernel=exponential_dispersal_kernel(distance_decay = 0.1),
                                                          stages = c(2,3)),
                                  mod = NULL,
                                  dens_dep = NULL)
  
  
  b_dynamics <- dynamics(habitat_dynamics = hab_dyn,
                         demography_dynamics = dem_dyn,
                         population_dynamics = pop_dyn)
  
  b_dynamics2 <- dynamics(habitat_dynamics = habitat_dynamics(),
                          demography_dynamics = demography_dynamics(),
                          population_dynamics = pop_dyn)
  
  b_dynamics3 <- dynamics(habitat_dynamics = habitat_dynamics(),
                          demography_dynamics = demography_dynamics(),
                          population_dynamics = pop_dyn2)
  
  b_dynamics4 <- dynamics(habitat_dynamics = habitat_dynamics(),
                          demography_dynamics = demography_dynamics(),
                          population_dynamics = pop_dyn3)
  
  b_dynamics5 <- dynamics(habitat_dynamics = hab_dyn,
                          demography_dynamics = dem_dyn,
                          population_dynamics = pop_dyn4)
  
  b_dynamics6 <- dynamics(habitat_dynamics = habitat_dynamics(),
                          demography_dynamics = dem_dyn2,
                          population_dynamics = pop_dyn)
  
  b_dynamics7 <- dynamics(habitat_dynamics = habitat_dynamics(),
                          demography_dynamics = dem_dyn3,
                          population_dynamics = pop_dyn)
  
  b_dynamics8 <- dynamics(habitat_dynamics = habitat_dynamics(),
                          demography_dynamics = demography_dynamics(),
                          population_dynamics = pop_dyn5)
  
  expect_true(inherits(simulation(state = b_state,
                                  dynamics = b_dynamics,
                                  timesteps = 5)[1],
                       "simulation_results"))
    
  expect_true(inherits(simulation(state = b_state,
                                  dynamics = b_dynamics,
                                  timesteps = 5),
                       "simulation_results"))
  
  expect_true(is.simulation_results(simulation(state = b_state,
                                  dynamics = b_dynamics,
                                  timesteps = 10)))

  expect_true(inherits(simulation(state = b_state2,
                                  dynamics = b_dynamics,
                                  timesteps = 10),
                       "simulation_results"))
    
  expect_true(inherits(simulation(state = b_state2,
                                  dynamics = b_dynamics2,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(state = b_state,
                                  dynamics = b_dynamics3,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(state = b_state2,
                                  dynamics = b_dynamics3,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(state = b_state,
                                  dynamics = b_dynamics4,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(state = b_state,
                                  dynamics = b_dynamics5,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(state = b_state,
                                  dynamics = b_dynamics6,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(state = b_state,
                                  dynamics = b_dynamics8,
                                  timesteps = 10),
                       "simulation_results"))

  expect_error(inherits(simulation(state = b_state,
                                   dynamics = b_dynamics7,
                                   timesteps = 10),
                        "simulation_results"))
  
  expect_error(simulation(state = b_state,
                                  dynamics = b_dynamics,
                                  timesteps = 15))
  
  expect_error(simulation(state = b_state3,
                          dynamics = b_dynamics,
                          timesteps = 10)
  )
  
  expect_error(simulation(state = b_state2,
                          dynamics = b_dynamics6,
                          timesteps = 10)
  )

  expect_error(simulation(state = b_state,
                          dynamics = b_dynamics,
                          timesteps = 15,
                          replicates = 5)
  )
  
  expect_error(simulation(state = b_state,
                          dynamics = b_dynamics6,
                          timesteps = 15,
                          replicates = 5)
  )
  
  expect_error(simulation(state = b_state3,
                          dynamics = b_dynamics,
                          timesteps = 10,
                          replicates = 5)
  )
  
  expect_true(inherits(simulation(state = b_state,
                                dynamics = b_dynamics,
                                timesteps = 10,
                                replicates = 2),
                       "simulation_results")
  )
   
  test_simulation <- simulation(state = b_state,
                                dynamics = b_dynamics,
                                timesteps = 10,
                                replicates = 5,
                                parallel = TRUE)
    
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
       object = "population",
       type = "graph",
       stage = 0)
  
  plot(test_simulation[1],
       object = "population",
       type = "graph",
       stage = 2)
  
  plot(test_simulation[1],
       object = "habitat_suitability")
  
  plot(test_simulation[1],
       object = "carrying_capacity")
    
  expect_error(plot(test_simulation,
                    object = "population",
                    type = "raster"))
  
  expect_error(plot(test_simulation,
                    object = "habitat_suitability"))
  
  expect_error(plot(test_simulation,
                    object = "carrying_capacity"))
  
})
