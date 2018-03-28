library(steps)
library(raster)
library(future)
library(rasterVis)

#system('rm /home/casey/Research/Github/steps/src/*.so /home/casey/Research/Github/steps/src/*.o')

#.pardefault <- par()
#par(mar=c(0.5,0.5,0.5,0.5))

############ MODEL INPUTS ###############

#Note, this is an age-structured matrix and is not generic
koala.trans.mat <- matrix(c(0.000,0.000,0.302,0.302,
                              0.940,0.000,0.000,0.000,
                              0.000,0.884,0.000,0.000,
                              0.000,0.000,0.793,0.793),
                            nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat) <- rownames(koala.trans.mat) <- c('Stage_0_1','Stage_1_2','Stage_2_3','Stage_3_')

koala.trans.mat.es <- matrix(c(0.000,0.000,1,1,
                                 1,0.000,0.000,0.000,
                                 0.000,1,0.000,0.000,
                                 0.000,0.000,1,1),
                               nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat.es) <- rownames(koala.trans.mat.es) <- c('Stage_0_1','Stage_1_2','Stage_2_3','Stage_3_')


koala.hab.suit <- raster("inst/extdata/Koala_HabSuit.tif") # read in spatial habitat suitability raster
koala.hab.suit <- (koala.hab.suit - cellStats(koala.hab.suit, min)) / (cellStats(koala.hab.suit, max) - cellStats(koala.hab.suit, min))
names(koala.hab.suit) <- "Habitat"
plot(koala.hab.suit, box = FALSE, axes = FALSE)

koala.hab.k <- koala.hab.suit*60
names(koala.hab.k) <- "Carrying Capacity"
plot(koala.hab.k, box = FALSE, axes = FALSE)

koala.pop <- stack(replicate(4, (koala.hab.k)*0.05))
names(koala.pop) <- colnames(koala.trans.mat)

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

koala.dist.fire <- stack(list.files("inst/extdata", full = TRUE, pattern = 'Koala_Fire*'))

####### Permutation 1 ########

koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k)
koala.demography <- build_demography(transition_matrix = koala.trans.mat,
                                     dispersal_parameters = rlnorm(1))
koala.population <- build_population(population_raster = koala.pop)
koala.state <- build_state(koala.habitat, koala.demography, koala.population)

koala.habitat.dynamics <- as.habitat_dynamics(no_habitat_dynamics)
koala.demography.dynamics <- as.demography_dynamics(no_demographic_dynamics)
koala.population.dynamics <- as.population_dynamics(fast_population_dynamics)
koala.dynamics <- build_dynamics(koala.habitat.dynamics,
                                 koala.demography.dynamics,
                                 koala.population.dynamics
)

######################################

####### Permutation 2 ########

koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k)
koala.demography <- build_demography(transition_matrix = koala.trans.mat,
                                     dispersal_parameters = koala.disp.param)
koala.population <- build_population(population_raster = koala.pop)
koala.state <- build_state(habitat = koala.habitat,
                           demography = koala.demography,
                           population = koala.population)

koala.habitat.dynamics <- as.habitat_dynamics(no_habitat_dynamics)
koala.demography.dynamics <- as.demography_dynamics(no_demographic_dynamics)
koala.population.dynamics <- as.population_dynamics(ca_dispersal_population_dynamics)
koala.dynamics <- build_dynamics(habitat_dynamics = koala.habitat.dynamics,
                                 demography_dynamics = koala.demography.dynamics,
                                 population_dynamics = koala.population.dynamics,
                                 order = c("habitat_dynamics",
                                           "demography_dynamics",
                                           "population_dynamics")
)


######################################

####### Permutation 3 ########

koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k)
koala.demography <- build_demography(transition_matrix = koala.trans.mat,
                                     dispersal_parameters = koala.disp.param)
koala.population <- build_population(population_raster = koala.pop)
koala.state <- build_state(habitat = koala.habitat,
                           demography = koala.demography,
                           population = koala.population)

koala.habitat.dynamics <- as.habitat_dynamics(no_habitat_dynamics)
koala.demography.dynamics <- envstoch_demography_dynamics(global_transition_matrix = koala.trans.mat,
                                                           stochasticity = koala.trans.mat.es)
koala.population.dynamics <- as.population_dynamics(ca_dispersal_population_dynamics)
koala.dynamics <- build_dynamics(habitat_dynamics = koala.habitat.dynamics,
                                 demography_dynamics = koala.demography.dynamics,
                                 population_dynamics = koala.population.dynamics,
                                 order = c("habitat_dynamics",
                                           "demography_dynamics",
                                           "population_dynamics")
)


######################################

####### Permutation 4 ########

koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k)
koala.demography <- build_demography(transition_matrix = koala.trans.mat,
                                     dispersal_parameters = koala.disp.param2)
koala.population <- build_population(population_raster = koala.pop)
koala.state <- build_state(habitat = koala.habitat,
                           demography = koala.demography,
                           population = koala.population)

koala.habitat.dynamics <- as.habitat_dynamics(no_habitat_dynamics)
koala.demography.dynamics <- envstoch_demography_dynamics(global_transition_matrix = koala.trans.mat,
                                                           stochasticity = koala.trans.mat.es)
koala.population.dynamics <- as.population_dynamics(ca_dispersal_population_dynamics)
koala.dynamics <- build_dynamics(habitat_dynamics = koala.habitat.dynamics,
                                 demography_dynamics = koala.demography.dynamics,
                                 population_dynamics = koala.population.dynamics,
                                 order = c("habitat_dynamics",
                                           "demography_dynamics",
                                           "population_dynamics")
)


######################################

####### Permutation 5 ########

koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k)
koala.demography <- build_demography(transition_matrix = koala.trans.mat,
                                     dispersal_parameters = koala.disp.param)
koala.population <- build_population(population_raster = koala.pop)
koala.state <- build_state(habitat = koala.habitat,
                           demography = koala.demography,
                           population = koala.population)

koala.habitat.dynamics <- fire_habitat_dynamics(habitat_suitability = koala.hab.suit,
                                                disturbance_layers = koala.dist.fire,
                                                effect_time=3)
koala.demography.dynamics <- envstoch_demography_dynamics(global_transition_matrix = koala.trans.mat,
                                                           stochasticity = koala.trans.mat.es)
koala.population.dynamics <- as.population_dynamics(ca_dispersal_population_dynamics)
koala.dynamics <- build_dynamics(habitat_dynamics = koala.habitat.dynamics,
                                 demography_dynamics = koala.demography.dynamics,
                                 population_dynamics = koala.population.dynamics,
                                 order = c("habitat_dynamics",
                                           "population_dynamics",
                                           "demography_dynamics")
)


#####################################

####### Permutation 6 ########

koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k*.07)
koala.demography <- build_demography(transition_matrix = koala.trans.mat,
                                     dispersal_parameters = koala.disp.param)
koala.population <- build_population(population_raster = koala.pop)
koala.state <- build_state(habitat = koala.habitat,
                           demography = koala.demography,
                           population = koala.population)

koala.habitat.dynamics <- fire_habitat_dynamics(habitat_suitability = koala.hab.suit,
                                                disturbance_layers = koala.dist.fire,
                                                effect_time=3)
koala.demography.dynamics <- envstoch_demography_dynamics(global_transition_matrix = koala.trans.mat,
                                                           stochasticity = koala.trans.mat.es)
koala.population.dynamics <- as.population_dynamics(fft_dispersal_population_dynamics)
koala.dynamics <- build_dynamics(habitat_dynamics = koala.habitat.dynamics,
                                 demography_dynamics = koala.demography.dynamics,
                                 population_dynamics = koala.population.dynamics,
                                 order = c("habitat_dynamics",
                                           "demography_dynamics",
                                           "population_dynamics")
)


#####################################

system.time(
  exp_results <- experiment(koala.state,
                         koala.dynamics,
                         timesteps = 20
                         )
)

plan(multiprocess) 
sim_results <- simulation(koala.state,
                          koala.dynamics,
                          timesteps = 20,
                          simulations = 5
                          )

plot(sim_results[[1]])

plot(exp_results)

plot(exp_results, stage = 2)

plot(exp_results, stage = 0)

plot(exp_results, type = "raster", stage = 2)

plot(exp_results, object = "habitat_suitability")

plot(exp_results, object = "carrying_capacity")

stages <- raster::nlayers(sim_results[[1]][[1]]$population$population_raster)

graph.pal <- c("#94d1c7",
               "#cccc2b",
               "#bebada",
               "#fb8072",
               "#80b1d3",
               "#fdb462",
               "#b3de69",
               "#fccde5",
               "#969696",
               "#bc80bd"
)

idx <- which(!is.na(raster::getValues(sim_results[[1]][[1]]$population$population_raster[[1]])))
pops <- lapply(sim_results, function(y) lapply(y, function(x) raster::extract(x$population$population_raster, idx)))
pop_sums <- lapply(pops, function(y) lapply(y, function(x) colSums(x)))

stage_names <- unlist(dimnames(sim_results[[1]][[1]]$demography$global_transition_matrix)[1])

graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,stages))

for (i in seq_len(stages)) {
  
  graphics::plot(unlist(lapply(pop_sums, function(y) lapply(y, function(x) x[[i]]))),
                 type='l',
                 ylab=paste("Total Population: ",stage_names[i]),
                 xlab="Time (years)",
                 lwd=2,
                 col=graph.pal[i]
  )
  graphics::abline(h=raster::cellStats(x[[1]]$habitat$carrying_capacity,sum)/stages,
                   lwd=1,
                   lty=2)
  
}