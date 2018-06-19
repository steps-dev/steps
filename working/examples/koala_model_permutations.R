library(raster)
library(future)
plan(multiprocess)

#system('rm /home/casey/Research/Github/steps/src/*.so /home/casey/Research/Github/steps/src/*.o')
#.pardefault <- par()
#par(mar=c(0.5,0.5,0.5,0.5))

############ MODEL INPUTS ###############

# Define a transition matrix. Note, this is an age-structured matrix and is not generic
koala.trans.mat <- matrix(c(0.000,0.000,0.302,0.302,
                            0.940,0.000,0.000,0.000,
                            0.000,0.884,0.000,0.000,
                            0.000,0.000,0.793,0.793),
                          nrow = 4, ncol = 4, byrow = TRUE)
# Add life-stage names to the matrix
colnames(koala.trans.mat) <- rownames(koala.trans.mat) <- c('Stage_0_1','Stage_1_2','Stage_2_3','Stage_3_')

# Define the standard deviation around the transition matrix values (for stochasticity).
# Note, can also be added as a single number in the function call.
koala.trans.mat.es <- matrix(c(0.000,0.000,1,1,
                               1,0.000,0.000,0.000,
                               0.000,1,0.000,0.000,
                               0.000,0.000,1,1),
                             nrow = 4, ncol = 4, byrow = TRUE)
# Add life-stage names to the stochasticity matrix
colnames(koala.trans.mat.es) <- rownames(koala.trans.mat.es) <- c('Stage_0_1','Stage_1_2','Stage_2_3','Stage_3_')

# Define an initial habitat suitability layer
koala.hab.suit <- raster("inst/extdata/Koala_HabSuit.tif") # read in spatial habitat suitability raster
koala.hab.suit <- (koala.hab.suit - cellStats(koala.hab.suit, min)) / # standardise values to be between 0 and 1
  (cellStats(koala.hab.suit, max) - cellStats(koala.hab.suit, min))
names(koala.hab.suit) <- "Habitat" # assign a name to the layer
plot(koala.hab.suit, box = FALSE, axes = FALSE)

# Define an initial carrying capacity layer
koala.hab.k <- ceiling(koala.hab.suit*20) # up to 20 individuals in the best habitat
names(koala.hab.k) <- "Carrying Capacity" # assign a name to the layer
plot(koala.hab.k, box = FALSE, axes = FALSE)

# Compute the stable age distributions for the total population
koala.stable.states <- abs(eigen(koala.trans.mat)$vectors[,1]/base::sum(eigen(koala.trans.mat)$vectors[,1]) ) 

# Define the initial populations using the stable age distributions and the carrying capacity
koala.pop <- ceiling(stack(replicate(4, koala.hab.k)) * koala.stable.states)
names(koala.pop) <- colnames(koala.trans.mat) # assign names to the layers
plot(koala.pop, box = FALSE, axes = FALSE)

# Define a landscape barrier to dispersal. Note, in this case a road.
koala.disp.bar <- koala.hab.suit*0
koala.disp.bar[cellFromRow(koala.disp.bar,nrow(koala.disp.bar)/2)] <- 1
plot(koala.disp.bar, box = FALSE, axes = FALSE)

# Define the initial dispersal parameters
koala.disp.param <- list(dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
                         dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
                         dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
                         barrier_type=1,
                         barriers_map=koala.disp.bar,
                         use_barriers=TRUE
)

# Load some habitat disturbance layers. Note, these must match the intended number of timesteps.
koala.dist.fire <- stack(list.files("inst/extdata", full = TRUE, pattern = 'Koala_Fire*'))

# Define some population disturbance layers. In this case where and how many individuals will be translocated.
koala.pop.source <- koala.pop[[4]]
koala.pop.source[] <- 0
koala.pop.source[sample(which(getValues(koala.pop[[4]]) >= 5), 5)] <- 1
plot(koala.pop.source, box = FALSE, axes = FALSE)

koala.pop.sink <- koala.pop[[4]]
koala.pop.sink[] <- 0
koala.pop.sink[sample(which(getValues(koala.pop[[1]]) == 1 |
                            getValues(koala.pop[[2]]) == 1 |
                            getValues(koala.pop[[3]]) == 1 |
                            getValues(koala.pop[[4]]) == 1),
                      cellStats(koala.pop.source, sum))] <- 1
plot(koala.pop.sink, box = FALSE, axes = FALSE)

# Note, for reintroductions from off-site, e.g. captive breeding, the source layer is all zeros and only the sink layer has values.
koala.pop.source.cb <- koala.pop.source
koala.pop.source.cb[] <- 0


# Define some demography disturbance layers. Note, these must match the intended number of timesteps.
koala.surv <- list(stack(list.files("inst/extdata", full = TRUE, pattern = 'Koala_Sur_F03R+')),
                   stack(list.files("inst/extdata", full = TRUE, pattern = 'Koala_Sur_F01+')),
                   stack(list.files("inst/extdata", full = TRUE, pattern = 'Koala_Sur_F02NR+')),
                   stack(list.files("inst/extdata", full = TRUE, pattern = 'Koala_Sur_F03NR+')))

koala.fec <- list(NULL,
                  NULL,
                  stack(list.files("inst/extdata", full = TRUE, pattern = 'Koala_Sur_F02R+')),
                  stack(list.files("inst/extdata", full = TRUE, pattern = 'Koala_Sur_F03R+')))

# Construct the habitat object
koala.habitat <- build_habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k,
                               misc = NULL)
# Construct the demography object
koala.demography <- build_demography(transition_matrix = koala.trans.mat,
                                     type = 'local', 
                                     habitat_suitability = koala.hab.suit,
                                     misc = NULL)
# Construct the population object
koala.population <- build_population(population_raster = koala.pop)

# Combine habitat, demography and population obects in a state object
koala.state <- build_state(habitat = koala.habitat,
                           demography = koala.demography,
                           population = koala.population)

# Create list of habitat dynamic modules
hab.dyn.list <- list(disturbance_fires(habitat_suitability = koala.hab.suit,
                                         disturbance_layers = koala.dist.fire,
                                         effect_time = 3),
                     NULL,
                     NULL)
# Create list of demography dynamic modules
dem.dyn.list <- list(demo_environmental_stochasticity(koala.trans.mat,
                                                      koala.trans.mat.es),
                     demo_density_dependence(koala.trans.mat,
                                             fecundity_fraction = 0.7,
                                             survival_fraction = 0.8),
                     surv_fec_modify(koala.trans.mat,
                                            koala.surv,
                                            koala.fec),
                     NULL)

# Create list of population dynamic modules
pop.dyn.list <- list(simple_growth(),
                     demographic_stochasticity(),
                     pop_density_dependence(),
                     simple_dispersal(dispersal_parameters = koala.disp.param),
                     kernel_function_dispersal(dispersal_parameters = koala.disp.param),
                     cellular_automata_dispersal(koala.disp.param),
                     fast_fourier_dispersal(koala.disp.param),
                     pop_translocation(koala.pop.source,
                                       koala.pop.sink,
                                       effect_timesteps = c(5,15)),
                     NULL)

# Create table with proposed combinations of dynamics for simulations
dyn_perm <- list(
  list(3, 4, 1, 5, 9, 9),
  list(3, 4, 1, 9, 9, 9),
  list(3, 4, 1, 3, 9, 9),
  list(3, 4, 1, 3, 7, 9),
  list(3, 4, 1, 3, 7, 9),
  list(3, 4, 2, 9, 9, 9),
  list(3, 4, 2, 3, 9, 9),
  list(3, 4, 2, 3, 7, 9),
  list(3, 4, 2, 3, 7, 9),
  list(3, 4, 2, 3, 7, 8),
  list(3, c(1,2), 1, 9, 9, 9),
  list(3, c(1,3), 1, 3, 9, 9),
  list(1, 3, 1, 3, 7, 9),
  list(3, 3, 1, 3, 7, 9),
  list(3, c(1,2), 2, 9, 9, 9),
  list(3, c(1,3), 2, 3, 9, 9),
  list(1, 3, 2, 3, 9, 9),
  list(3, 3, 2, 3, 7, 9),
  list(1, 1, 2, 3, 7, 8)
)

for ( i in seq_len(length(dyn_perm)) ) {
  koala.state <- build_state(habitat = koala.habitat,
                             demography = koala.demography,
                             population = koala.population)
  
  
  koala.dynamics <- build_dynamics(habitat_dynamics = do.call(habitat_dynamics, hab.dyn.list[dyn_perm[[i]][[1]]]),
                                   demography_dynamics = do.call(demography_dynamics, dem.dyn.list[dyn_perm[[i]][[2]]]),
                                   population_dynamics = population_dynamics(pop_change = pop.dyn.list[[dyn_perm[[i]][[3]]]],
                                                                             pop_dens_dep = pop.dyn.list[[dyn_perm[[i]][[4]]]],
                                                                             pop_disp = pop.dyn.list[[dyn_perm[[i]][[5]]]],
                                                                             pop_mod = pop.dyn.list[[dyn_perm[[i]][[6]]]]),
                                   order = c("habitat_dynamics",
                                             "population_dynamics",
                                             "demography_dynamics")
  )
  
  
  sim_results <- simulation(state = koala.state,
                            dynamics = koala.dynamics,
                            timesteps = 10,
                            replicates = 1)
  
  pdf(paste0("working/examples/plots/pop_perm_",sprintf("%03d", i),".pdf"))
  plot(sim_results)
  title(main = paste0("Model Permutation ",sprintf("%03d", i)))
  dev.off()
}


plot(sim_results)

plot(sim_results[[1]][[10]]$population$population_raster)
