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
cells <- sample(c(1:ncell(r2)), 10)
r2[c(adjacent(hab.suit, cells, directions=8, pairs=FALSE),cells)]  <- 10
r3 <- r2#*hab.suit

pop <- stack(r3*1,r3*2,r3*3,r3*2)

hab.k <- hab.suit*10

params <- list(
  dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=1,'Stage_2-3'=1,'Stage_3+'=0),
  dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
  dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0))

b_hab <- build_habitat(habitat_suitability = hab.suit,
                       carrying_capacity = hab.k)
b_pop <- build_population(pop)
b_dem <- build_demography(transition_matrix = mat,
                          dispersal_parameters = params)

b_state <- build_state(habitat = b_hab,
                       population = b_pop,
                       demography = b_dem)

hab_dyn <- no_habitat_dynamics()

dem_dyn <- demography_dynamics(env_stoch = demo_environmental_stochasticity(global_transition_matrix = mat,
                                                                        stochasticity = mat_sd),
                               demo_dens_dep = demo_density_dependence(fecundity_fraction = 1,
                                                                   survival_fraction = 0.5))
dem_dyn2 <- demography_dynamics()

pop_dyn <- ca_dispersal_population_dynamics()

b_dynamics <- build_dynamics(habitat_dynamics = hab_dyn,
                             demography_dynamics = dem_dyn2,
                             population_dynamics = pop_dyn)

result <- simulation(state = b_state,
           dynamics = b_dynamics,
           timesteps = 10,
           replicates = 10)

plot(result)
