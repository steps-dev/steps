context('simulation_results-class')

test_that('simulation_results classes work', {
  library(raster)
  library(rgdal)
  library(future)
  plan(multiprocess)

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
  cells <- sample(c(1:ncell(r2)), 10)
  r2[c(adjacent(hab.suit, cells, directions=8, pairs=FALSE),cells)]  <- 10
  r3 <- r2#*hab.suit
  
  pop <- stack(r3*2,r3*2,r3*2,r3*2)
  
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
  b_pop <- build_population(pop)
  b_dem <- build_demography(transition_matrix = mat,
                          dispersal_parameters = params2)
  
  b_hab2 <- build_habitat(habitat_suitability = hab.suit,
                         carrying_capacity = NULL)
  
  b_dem2 <- build_demography(transition_matrix = fake_transition_matrix(4),
                             dispersal_parameters = params3)
  
  b_dem3 <- build_demography(transition_matrix = mat,
                            dispersal_parameters = rlnorm(1))
  
  b_state <- build_state(habitat = b_hab,
              population = b_pop,
              demography = b_dem)
  
  b_state2 <- build_state(habitat = b_hab,
                         population = b_pop,
                         demography = b_dem2)
  
  b_state3 <- build_state(habitat = b_hab2,
                          population = b_pop,
                          demography = b_dem)
  
  b_state4 <- build_state(habitat = b_hab,
                          population = b_pop,
                          demography = b_dem3)
  
  hab_dyn <- fire_habitat_dynamics(habitat_suitability = hab.suit,
                                                  disturbance_layers = dist.s,
                                                  effect_time = 2)
  dem_dyn <- envstoch_demography_dynamics(global_transition_matrix = mat,
                                           stochasticity = mat_sd)
  pop_dyn <- ca_dispersal_population_dynamics()
  pop_dyn2 <- fast_population_dynamics()
  pop_dyn3 <- fft_dispersal_population_dynamics()
  
  b_dynamics <- build_dynamics(habitat_dynamics = hab_dyn,
                               demography_dynamics = dem_dyn,
                               population_dynamics = pop_dyn)
  
  b_dynamics2 <- build_dynamics(habitat_dynamics = no_habitat_dynamics(),
                               demography_dynamics = no_demography_dynamics(),
                               population_dynamics = pop_dyn)
  
  b_dynamics3 <- build_dynamics(habitat_dynamics = no_habitat_dynamics(),
                                demography_dynamics = no_demography_dynamics(),
                                population_dynamics = pop_dyn2)
  
  b_dynamics4 <- build_dynamics(habitat_dynamics = hab_dyn,
                               demography_dynamics = dem_dyn,
                               population_dynamics = pop_dyn3)

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
                                  dynamics = b_dynamics2,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(state = b_state4,
                                  dynamics = b_dynamics3,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_true(inherits(simulation(state = b_state,
                                  dynamics = b_dynamics4,
                                  timesteps = 10),
                       "simulation_results"))
  
  expect_error(simulation(state = b_state,
                                  dynamics = b_dynamics,
                                  timesteps = 15))
  
  expect_error(simulation(state = b_state3,
                          dynamics = b_dynamics,
                          timesteps = 10)
  )

  expect_error(simulation(state = b_state,
                          dynamics = b_dynamics,
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