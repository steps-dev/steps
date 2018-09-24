## ---- message = FALSE, echo = FALSE--------------------------------------
library(steps)
library(raster)

## ---- message = FALSE----------------------------------------------------

koala.trans.mat <- matrix(c(0.000,0.000,0.302,0.302,
                              0.940,0.000,0.000,0.000,
                              0.000,0.884,0.000,0.000,
                              0.000,0.000,0.793,0.793),
                            nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat) <- rownames(koala.trans.mat) <- c('Juveniles','Sub_Adults','Adults','Super_Adults')

koala.trans.mat.es <- matrix(c(0.000,0.000,1,1,
                                 1,0.000,0.000,0.000,
                                 0.000,1,0.000,0.000,
                                 0.000,0.000,1,1),
                               nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat.es) <- rownames(koala.trans.mat.es) <- c('Juveniles','Sub_Adults','Adults','Super_Adults')


## ---- message = FALSE, fig.align="center"--------------------------------

# read in spatial habitat suitability raster
koala.hab.suit <- raster(system.file("extdata","Koala_HabSuit.tif", package="steps"))

# rescale habitat raster to relative suitability (or likelihood of occurrence) between 0 and 1.
koala.hab.suit <- (koala.hab.suit - cellStats(koala.hab.suit, min)) / (cellStats(koala.hab.suit, max) - cellStats(koala.hab.suit, min))

# rename raster so that it is tracked throughout the simulations
names(koala.hab.suit) <- "Habitat"

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(koala.hab.suit, box = FALSE, axes = FALSE)

## ---- message = FALSE, fig.align="center"--------------------------------

# create carrying capacity layer using the habitat suitability raster (or provide a custom one)
koala.hab.k <- ceiling(koala.hab.suit * 60)

# rename raster so that it is tracked throughout the simulations
names(koala.hab.k) <- "Carrying Capacity"

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(koala.hab.k, box = FALSE, axes = FALSE)

## ---- message = FALSE, fig.align="center"--------------------------------

# create population layers using the carrying capacity raster (or provide a custom ones)
koala.pop <- stack(replicate(4, ceiling(koala.hab.k * 0.2)))

# rename stack so that the layers are tracked throughout the simulations
names(koala.pop) <- colnames(koala.trans.mat)

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spplot(koala.pop)

## ---- message = FALSE----------------------------------------------------

# define all of the dispersal distances
dispersal_distance <- list('Juveniles'=0,'Sub_Adults'=10,'Adults'=10,'Super_Adults'=0)

# define the dispersal kernals for each life-stage
dispersal_kernel <- list('Juveniles'=0,'Sub_Adults'=exp(-c(0:9)^1/3.36),'Adults'=exp(-c(0:9)^1/3.36),'Super_Adults'=0)

# define all of the dispersal proportions
dispersal_proportion <- list('Juveniles'=0,'Sub_Adults'=0.35,'Adults'=0.35*0.714,'Super_Adults'=0)

# combine all of the parameters in a list
koala.disp.param <- list(dispersal_distance=dispersal_distance,
                                      dispersal_kernel=dispersal_kernel,
                                      dispersal_proportion=dispersal_proportion
                         )

## ---- message = FALSE----------------------------------------------------
koala.dist.fire <- stack(list.files(system.file("extdata", package="steps"), full.names = TRUE, pattern = 'Koala_Fire*'))

## ---- message = FALSE----------------------------------------------------
koala.habitat <- habitat(habitat_suitability = koala.hab.suit,
                               carrying_capacity = koala.hab.k)

## ---- message = FALSE----------------------------------------------------
koala.demography <- demography(transition_matrix = koala.trans.mat,
                                     type = 'local',
                                     habitat_suitability = koala.hab.suit)

## ---- message = FALSE----------------------------------------------------
koala.population <- population(population_raster = koala.pop)

## ---- message = FALSE----------------------------------------------------
koala.state <- state(habitat = koala.habitat,
                           demography = koala.demography,
                           population = koala.population)

## ---- message = FALSE----------------------------------------------------
koala.habitat.dynamics <- habitat_dynamics(disturbance_fires(habitat_suitability = koala.hab.suit,
                                                disturbance_layers = koala.dist.fire,
                                                effect_time=3))

## ---- message = FALSE----------------------------------------------------
koala.demography.dynamics <- demography_dynamics(environmental_stochasticity(transition_matrix = koala.trans.mat,
                                                                                  stochasticity = koala.trans.mat.es),
                                                 density_dependence(transition_matrix = koala.trans.mat))


## ---- message = FALSE----------------------------------------------------
koala.population.dynamics <- population_dynamics(change = simple_growth(),
                                                 disp = cellular_automata_dispersal(dispersal_distance = list(0,
                                                                                                           10,
                                                                                                           10,
                                                                                                           0),
                                                                                        dispersal_kernel = list(0,
                                                                                                         exp(-c(0:9)^1/3.36),
                                                                                                         exp(-c(0:9)^1/3.36),
                                                                                                         0),
                                                                                        dispersal_proportion = list(0,
                                                                                                             0.35,
                                                                                                             0.35*0.714,
                                                                                                             0),
                                                                                        barrier_type = 0,
                                                                                        dispersal_steps = 1,
                                                                                        use_barriers = FALSE,
                                                                                        barriers_map = NULL,
                                                                                        arrival_probability = "habitat_suitability",
                                                                                        carrying_capacity = "carrying_capacity"),
                                          mod = NULL,
                                          dens_dep = ceiling_density_dependence(stages = c(3,4)))

## ---- message = FALSE----------------------------------------------------
koala.dynamics <- dynamics(habitat_dynamics = koala.habitat.dynamics,
                                 demography_dynamics = koala.demography.dynamics,
                                 population_dynamics = koala.population.dynamics,
                                 order = c("habitat_dynamics",
                                           "demography_dynamics",
                                           "population_dynamics")
)

## ---- message = FALSE,  results='hide'-----------------------------------
koala.results <- simulation(koala.state,
                       koala.dynamics,
                       timesteps = 20,
                       replicates = 1
                       )

## ---- message = FALSE, fig.width=7, fig.align="center"-------------------
plot(koala.results)

## ---- message = FALSE, fig.width=4, fig.align="center"-------------------
plot(koala.results, stage = 0)

## ---- message = FALSE, fig.width=4, fig.align="center"-------------------
plot(koala.results, stage = 2, newplot = TRUE)

## ---- message = FALSE, fig.width=7, fig.align="center"-------------------
plot(koala.results, type = "raster", stage = 2)

## ---- message = FALSE, fig.width=7, fig.align="center"-------------------
plot(koala.results, object = "habitat_suitability")

plot(koala.results, object = "carrying_capacity")

## ---- message = FALSE, echo = FALSE--------------------------------------
koala.sim.results <- simulation(koala.state,
                                koala.dynamics,
                                timesteps = 20,
                                replicates = 3,
                                parallel = TRUE
                                )

## ---- message = FALSE, fig.width=4, fig.align="center"-------------------
plot(koala.sim.results)

## ---- message = FALSE, fig.width=4, fig.align="center"-------------------
plot(koala.sim.results[3], stage = 0)

## ---- message = FALSE, fig.width=7, fig.align="center"-------------------
plot(koala.sim.results[1], type = "raster", stage = 2)

