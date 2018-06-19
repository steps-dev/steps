context('simulation_results-class')

test_that('simulation_results classes work', {
  
  library(raster)
  library(rgdal)
  library(future)
  plan(sequential)

  # the types of demography
  mat <- matrix(c(0.000,0.000,0.302,0.302,
                  0.940,0.000,0.000,0.000,
                  0.000,0.884,0.000,0.000,
                  0.000,0.000,0.793,0.793),
                nrow = 4, ncol = 4, byrow = TRUE)
  colnames(mat) <- rownames(mat) <- c('Stage_0-1','Stage_1-2','Stage_2-3','Stage_3+')
  
  # the types of demography
  mat_sd <- matrix(c(0.000,0.000,1,1,
                     1,0.000,0.000,0.000,
                     0.000,1,0.000,0.000,
                     0.000,0.000,1,1),
                   nrow = 4, ncol = 4, byrow = TRUE)
  colnames(mat_sd) <- rownames(mat_sd) <- c('Stage_0-1','Stage_1-2','Stage_2-3','Stage_3+')

  # the types of habitat attributes
  r <- raster(vals=1, nrows=150, ncols=150, res=c(5,5), crs=('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  hab.suit <- r*sample(seq(0,1,.01), ncell(r), replace=TRUE)
  
  r2 <- r
  r2[] <- 0
  cells <- sample(c(1:ncell(r2)), 150)
  r2[c(adjacent(hab.suit, cells, directions=16, pairs=FALSE),cells)]  <- 50
  r3 <- r2*hab.suit
  
  pop <- stack(r3*2,r3*3,r3*2,r3*3)
  
  hab.k <- hab.suit*10
  
  disp.bar <- hab.suit*0
  disp.bar[cellFromCol(disp.bar,ncol(disp.bar)/2)] <- 1
  disp.bar2 <- hab.suit*0
  disp.bar2[sampleRandom(disp.bar2, size=700, na.rm=TRUE, sp=TRUE)] <- 1
  
  dist_list <- list()
  
  for (i in 1:10) {
    
    r2 <- r
    r2[] <- 1
    cells <- sample(c(1:ncell(r2)), 10)
    r2[c(adjacent(hab.suit, cells, directions=8, pairs=FALSE),cells)]  <- 0.5
    dist_list[[i]] <- r2
    
  }
  
  dist.s <- stack(dist_list)

  surv_fec <- list(dist.s, dist.s, dist.s, dist.s)
  surv_fec2 <- list(dist.s, dist.s, dist.s, NULL)
  
  pop_source <- pop[[3]]
  pop_source[] <- 0
  pop_source[sample(which(getValues(pop[[3]]) >= 50), 3)] <- 1
  #plot(pop_source, box = FALSE, axes = FALSE)
  
  pop_sink <- pop[[3]]
  pop_sink[] <- 0
  pop_sink[sample(which(getValues(pop[[1]]) != 0 |
                        getValues(pop[[2]]) != 0 |
                        getValues(pop[[3]]) != 0 |
                        getValues(pop[[4]]) == 0),
                        cellStats(pop_source, sum))] <- 1
  #plot(pop_sink, box = FALSE, axes = FALSE)

  params <- list(
    dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=1,'Stage_2-3'=1,'Stage_3+'=0),
    dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
    dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0))
  
  params2 <- list(
    dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
    dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
    dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
    barriers_map=disp.bar,
    use_barriers=TRUE,
    barrier_type=1)
  
  params3 <- list(
    dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
    dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
    dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
    barriers_map=disp.bar2,
    use_barriers=TRUE,
    barrier_type=1)
  
  b_hab <- build_habitat(habitat_suitability = hab.suit,
                         carrying_capacity = hab.k)
  
  b_hab2 <- build_habitat(habitat_suitability = hab.suit,
                          carrying_capacity = NULL)
  
  b_pop <- build_population(pop)
 
  b_dem <- build_demography(transition_matrix = mat,
                            type = "local",
                            habitat_suitability = hab.suit)
 
  b_dem2 <- build_demography(transition_matrix = mat)
  
  b_state <- build_state(habitat = b_hab,
                         population = b_pop,
                         demography = b_dem)
  
  b_state2 <- build_state(habitat = b_hab,
                          population = b_pop,
                          demography = b_dem2)
  
  b_state3 <- build_state(habitat = b_hab2,
                          population = b_pop,
                          demography = b_dem)
  
  hab_dyn <- habitat_dynamics(disturbance_fires(habitat_suitability = hab.suit,
                                                  disturbance_layers = dist.s,
                                                  effect_time = 2))
  
  dem_dyn <- demography_dynamics(demo_environmental_stochasticity(transition_matrix = mat,
                                                                  stochasticity = mat_sd),
                                 demo_density_dependence(transition_matrix = mat,
                                                         fecundity_fraction = 0.8,
                                                         survival_fraction = 0.8))
    
  dem_dyn2 <- demography_dynamics(surv_fec_modify(transition_matrix = mat,
                                                  surv_layers = surv_fec,
                                                  fec_layers = surv_fec))
  
  dem_dyn3 <- demography_dynamics(surv_fec_modify(transition_matrix = mat,
                                                  surv_layers = surv_fec2,
                                                  fec_layers = surv_fec))

  pop_dyn <- population_dynamics(pop_change = simple_growth(),
                                 pop_disp = cellular_automata_dispersal(params2),
                                 pop_mod = pop_translocation(source_layer = pop_source,
                                                             sink_layer = pop_sink,
                                                             stages = NULL,
                                                             effect_timesteps = 2),
                                 pop_dens_dep = pop_density_dependence())
  pop_dyn2 <- population_dynamics(pop_change = demographic_stochasticity(),
                                  pop_disp = simple_dispersal(dispersal_parameters = params),
                                  pop_mod = NULL,
                                  pop_dens_dep = NULL)
  pop_dyn3 <- population_dynamics(pop_change = demographic_stochasticity(),
                                  pop_disp = fast_fourier_dispersal(params),
                                  pop_mod = NULL,
                                  pop_dens_dep = NULL)
  pop_dyn4 <- population_dynamics(pop_change = simple_growth(),
                                  pop_disp = cellular_automata_dispersal(params3),
                                  pop_mod = pop_translocation(source_layer = pop_source,
                                                             sink_layer = pop_sink,
                                                             stages = 4,
                                                             effect_timesteps = 2),
                                  pop_dens_dep = pop_density_dependence())

  b_dynamics <- build_dynamics(habitat_dynamics = hab_dyn,
                               demography_dynamics = dem_dyn,
                               population_dynamics = pop_dyn)
  
  b_dynamics2 <- build_dynamics(habitat_dynamics = habitat_dynamics(),
                                demography_dynamics = demography_dynamics(),
                                population_dynamics = pop_dyn)
  
  b_dynamics3 <- build_dynamics(habitat_dynamics = habitat_dynamics(),
                                demography_dynamics = demography_dynamics(),
                                population_dynamics = pop_dyn2)
  
  b_dynamics4 <- build_dynamics(habitat_dynamics = habitat_dynamics(),
                                demography_dynamics = demography_dynamics(),
                                population_dynamics = pop_dyn3)
  
  b_dynamics5 <- build_dynamics(habitat_dynamics = hab_dyn,
                                demography_dynamics = dem_dyn,
                                population_dynamics = pop_dyn4)
  
  b_dynamics6 <- build_dynamics(habitat_dynamics = habitat_dynamics(),
                                demography_dynamics = dem_dyn2,
                                population_dynamics = pop_dyn)
  
  b_dynamics7 <- build_dynamics(habitat_dynamics = habitat_dynamics(),
                                demography_dynamics = dem_dyn3,
                                population_dynamics = pop_dyn)
  
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
                                replicates = 5),
                       "simulation_results")
  )
   
  test_simulation <- simulation(state = b_state,
                                dynamics = b_dynamics,
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
