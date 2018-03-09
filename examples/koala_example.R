#library(dhmpr)
library(raster)

#system('rm /home/casey/Research/Github/dhmpr/src/*.so /home/casey/Research/Github/dhmpr/src/*.o')

#.pardefault <- par()
#par(mar=c(0.5,0.5,0.5,0.5))

############ MODEL INPUTS ###############

#Note, this is an age-structured matrix and is not generic
koala.trans.mat <- matrix(c(0.000,0.000,0.302,0.302,
                              0.940,0.000,0.000,0.000,
                              0.000,0.884,0.000,0.000,
                              0.000,0.000,0.793,0.793),
                            nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat) <- rownames(koala.trans.mat) <- c('Stage_0-1','Stage_1-2','Stage_2-3','Stage_3+')

koala.trans.mat.es <- matrix(c(0.000,0.000,1,1,
                                 1,0.000,0.000,0.000,
                                 0.000,1,0.000,0.000,
                                 0.000,0.000,1,1),
                               nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat.es) <- rownames(koala.trans.mat.es) <- c('Stage_0-1','Stage_1-2','Stage_2-3','Stage_3+')


koala.hab.suit <- raster("inst/extdata/koala/habitat/HS_crop_aggregate.tif") # read in spatial habitat suitability raster
koala.hab.suit <- (koala.hab.suit - cellStats(koala.hab.suit, min)) / (cellStats(koala.hab.suit, max) - cellStats(koala.hab.suit, min))
plot(koala.hab.suit, box = FALSE, axes = FALSE)

koala.hab.k <- koala.hab.suit*60
plot(koala.hab.k, box = FALSE, axes = FALSE)

koala.pop <- stack(replicate(4, (koala.hab.k)*0.05))

koala.disp.param <- list(dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
                                      dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
                                      dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0)
                         )

koala.disp.bar <- koala.hab.suit*0
koala.disp.bar[cellFromRow(koala.disp.bar,nrow(koala.disp.bar)/2)] <- 1
plot(koala.disp.bar, box = FALSE, axes = FALSE)

koala.disp.param2 <- list(dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
                         dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
                         dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
                         barrier_type=1,
                         barriers_map=koala.disp.bar,
                         use_barriers=TRUE
)

koala.disp.bar2 <- koala.hab.suit*0
koala.disp.bar2[sampleRandom(koala.disp.bar2, size=1000, na.rm=TRUE, sp=TRUE)] <- 1

koala.disp.param3 <- list(dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
                          dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
                          dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
                          barrier_type=1,
                          barriers_map=koala.disp.bar2,
                          use_barriers=TRUE
)

####### Permutation 1 ########

koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k)
koala.demography <- build_demography(koala.trans.mat,
                                     rlnorm(1))
koala.population <- build_population(population_raster = koala.pop)
koala.state <- build_state(koala.habitat, koala.demography, koala.population)

koala.habitat.dynamics <- as.habitat_dynamics(no_habitat_dynamics)
koala.demography.dynamics <- as.demography_dynamics(no_demographic_dynamics)
koala.population.dynamics <- as.population_dynamics(fast_population_dynamics)
koala.dynamics <- build_dynamics(koala.habitat.dynamics,
                                 koala.demography.dynamics,
                                 koala.population.dynamics
)

my.results <- experiment(koala.state,
                         koala.dynamics,
                         timesteps = 5)

rasters <- lapply(my.results, function (state) state$population$population_raster[[1]])
plot(stack(rasters))

######################################

####### Permutation 2 ########

koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k)
koala.demography <- build_demography(koala.trans.mat,
                                     koala.disp.param3)
koala.population <- build_population(population_raster = koala.pop)
koala.state <- build_state(koala.habitat, koala.demography, koala.population)

koala.habitat.dynamics <- as.habitat_dynamics(no_habitat_dynamics)
koala.demography.dynamics <- as.demography_dynamics(no_demographic_dynamics)
koala.population.dynamics <- as.population_dynamics(ca_dispersal_population_dynamics)
koala.dynamics <- build_dynamics(koala.habitat.dynamics,
                                 koala.demography.dynamics,
                                 koala.population.dynamics
)

my.results <- experiment(koala.state,
                         koala.dynamics,
                         timesteps = 10)

rasters <- lapply(my.results, function (state) state$population$population_raster[[2]])
plot(stack(rasters))

######################################

####### Permutation 3 ########

koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k)
koala.demography <- build_demography(koala.trans.mat,
                                     koala.disp.param)
koala.population <- build_population(population_raster = koala.pop)
koala.state <- build_state(koala.habitat, koala.demography, koala.population)

koala.habitat.dynamics <- as.habitat_dynamics(no_habitat_dynamics)
koala.demography.dynamics <- as.demography_dynamics(no_demographic_dynamics)
koala.population.dynamics <- as.population_dynamics(fft_dispersal_population_dynamics)
koala.dynamics <- build_dynamics(koala.habitat.dynamics,
                                 koala.demography.dynamics,
                                 koala.population.dynamics
)

my.results <- experiment(koala.state,
                         koala.dynamics,
                         timesteps = 10)

rasters <- lapply(my.results, function (state) state$population$population_raster[[2]])
plot(stack(rasters))

######################################