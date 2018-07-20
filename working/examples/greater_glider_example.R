library(steps)
library(raster)

############ MODEL INPUTS ###############

#Note, this is an age-structured matrix and is not generic
gg.trans.mat <- matrix(c(0.00,0.00,0.50,
                         0.50,0.00,0.00,
                         0.00,0.85,0.85),
                         nrow = 3, ncol = 3, byrow = TRUE)
colnames(gg.trans.mat) <- rownames(gg.trans.mat) <- c('Newborn','Juvenile','Adult')

gg.stable.states <- abs(eigen(gg.trans.mat)$vectors[,1]/base::sum(eigen(gg.trans.mat)$vectors[,1]) ) 

#r <- raster(vals=1, ymx=5860000, ymn=5800000, xmx=410000, xmn=360000, res=c(125,125), crs=('+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
r <- raster(vals=1, ymx=5860000, ymn=5800000, xmx=410000, xmn=360000, res=c(500,500), crs=('+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))

# tree <- resample(crop(raster("inst/extdata/Tree_Density_150.tif")/100, r),r)
# elev <- resample(crop(raster("inst/extdata/Elev_125.tif"), r),r)
# elev <- ((elev-cellStats(elev,"min"))/(cellStats(elev,"max")-cellStats(elev,"min")))
# gg.hab.suit <- r * tree * elev

gg.hab.suit <- resample(crop(raster("inst/extdata/GG_habsuit_125.tif"), r),r)
gg.hab.suit <- ((gg.hab.suit-cellStats(gg.hab.suit,"min"))/(cellStats(gg.hab.suit,"max")-cellStats(gg.hab.suit,"min")))
names(gg.hab.suit) <- "Habitat"
#plot(gg.hab.suit, box = FALSE, axes = FALSE, col = viridis(100))

#gg.hab.k <- ceiling(gg.hab.suit*3)
gg.hab.k <- ceiling(gg.hab.suit*48)
names(gg.hab.k) <- "Carrying Capacity"
#plot(gg.hab.k, box = FALSE, axes = FALSE, col = viridis(100))

# gg.pop <- floor(stack(replicate(3, gg.hab.k)) * gg.stable.states)

# gg.ala.occurrences <- read.csv("inst/extdata/ALA_GG_2018-06-26.csv")[,c(58:60,92:93,95)]
# coordinates(gg.ala.occurrences) <- ~decimalLongitude+decimalLatitude
# proj4string(gg.ala.occurrences) <- CRS("+proj=longlat +datum=WGS84")
# gg.ala.occurrences.r <- rasterize(crop(spTransform(gg.ala.occurrences, CRS("+init=epsg:28355")), r), r, field=1, fun=max)
# gg.pop <- ceiling(stack(replicate(3, gg.ala.occurrences.r)) * gg.stable.states)

#gg.popN <- stack(replicate(3, (gg.hab.suit*3)))
gg.popN <- stack(replicate(ncol(gg.trans.mat), (gg.hab.suit*48)))
gg.popN <-gg.popN*gg.stable.states

gg.pop<-stack()
for(i in 1:nlayers(gg.popN)){
  m <- ceiling(cellStats(gg.popN[[i]], max, na.rm=T))
  gg.pop1 <- calc(gg.popN[[i]], fun=function(x) rbinom(prob =(x/m), size=m, n=1))
  gg.pop <- stack(gg.pop,gg.pop1)
}
TotpopN <- sum(cellStats(gg.pop, 'sum', na.rm=T)) # Get total population size to check sensible
init.pop.size <- sum(gg.pop)
#plot(init.pop.size, box = FALSE, axes = FALSE, col = viridis(100))

names(gg.pop) <- colnames(gg.trans.mat)
#plot(gg.pop, box = FALSE, axes = FALSE, col = viridis(100))

gg.disp.bar <- gg.hab.suit*0
gg.disp.bar[cellFromRow(gg.disp.bar,c(nrow(gg.disp.bar)/2,(nrow(gg.disp.bar)/2)+1))] <- 1
#plot(gg.disp.bar, box = FALSE, axes = FALSE, col = viridis(100))

# koala.dist.fire <- stack(list.files("inst/extdata", full = TRUE, pattern = 'Koala_Fire*'))
# 
# gg.dist.fire <- stack()
# for (i in 1:nlayers(koala.dist.fire)) {
#   gg.fire <- r
#   gg.fire[] <- koala.dist.fire[[i]][]
#   gg.dist.fire <- stack(gg.dist.fire, gg.fire)
# }
# 
# gg.dist.fire2 <- stack(replicate(5, unlist(gg.dist.fire)))

# 
# gg.pop.source <- gg.pop[[4]]
# gg.pop.source[] <- 0
# gg.pop.source[sample(which(getValues(gg.pop[[4]]) >= 3), 25)] <- 1
# plot(gg.pop.source, box = FALSE, axes = FALSE)
# 
# gg.pop.sink <- gg.pop[[4]]
# gg.pop.sink[] <- 0
# gg.pop.sink[sample(which(getValues(gg.pop[[1]]) == 1 |
#                             getValues(gg.pop[[2]]) == 1 |
#                             getValues(gg.pop[[3]]) == 1 |
#                             getValues(gg.pop[[4]]) == 1),
#                       cellStats(gg.pop.source, sum))] <- 1
# plot(gg.pop.sink, box = FALSE, axes = FALSE)
# 
# gg.surv <- list(stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F03R+')),
#                    stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F01+')),
#                    stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F02NR+')),
#                    stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F03NR+')))
# 
# gg.fec <- list(NULL,
#                   NULL,
#                   stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F02R+')),
#                   stack(list.files("inst/extdata", full = TRUE, pattern = 'gg_Sur_F03R+')))



gg.habitat <- build_habitat(habitat_suitability = gg.hab.suit,
                            carrying_capacity = gg.hab.k)
gg.demography <- build_demography(transition_matrix = gg.trans.mat,
                                  scale = 'local',
                                  habitat_suitability = gg.hab.suit
                                  )
gg.population <- build_population(population_raster = gg.pop)
gg.state <- build_state(habitat = gg.habitat,
                        demography = gg.demography,
                        population = gg.population)


################# DETERMINISTIC GROWTH WITH ENV/DEMO STOCHASTICITY AND CA DISPERSAL #####################

gg.habitat.dynamics <- build_habitat_dynamics(#disturbance_fires(habitat_suitability = gg.hab.suit,
                                                         # disturbance_layers = gg.dist.fire2,
                                                          #effect_time=5)
                                        )
gg.demography.dynamics <- build_demography_dynamics(demo_environmental_stochasticity(transition_matrix = gg.trans.mat,
                                                                               stochasticity = 0.1)
                                              )
gg.population.dynamics <- build_population_dynamics(pop_change = simple_growth(demo_stoch = FALSE),
                                              #pop_disp = cellular_automata_dispersal(dispersal_distance=list(0, 16, 0),
                                                                          #dispersal_kernel=list(0, exp(-c(0:19)^1/10), 0),
                                                                          #dispersal_proportion=list(0, 0.5, 0)),
                                              pop_disp = kernel_dispersal(dispersal_proportion = list(0, 0.5, 0),
                                                                          dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.5),
                                                                          arrival_probability = NULL,
                                                                          fft = TRUE),
                                              pop_dens_dep = pop_density_dependence(stages = c(2,3))
                                              )
gg.dynamics <- build_dynamics(habitat_dynamics = gg.habitat.dynamics,
                              demography_dynamics = gg.demography.dynamics,
                              population_dynamics = gg.population.dynamics,
                              order = c("habitat_dynamics",
                                        "population_dynamics",
                                        "demography_dynamics")
)

sim_results <- simulation(state = gg.state,
                          dynamics = gg.dynamics,
                          timesteps = 100,
                          replicates = 3,
                          parallel = TRUE)

plot(sim_results)
########################################################################

# ################# WITH DEMOGRAPHIC DENSITY DEPENDENCE #####################
# 
# gg.habitat <- build_habitat(habitat_suitability = gg.hab.suit,
#                                carrying_capacity = gg.hab.k*0.1,
#                                misc = NULL)
# gg.demography <- build_demography(transition_matrix = gg.trans.mat,
#                                      scale = 'local', 
#                                      habitat_suitability = gg.hab.suit,
#                                      dispersal_parameters = gg.disp.param,
#                                      misc = NULL)
# gg.population <- build_population(population_raster = gg.pop)
# gg.state <- build_state(habitat = gg.habitat,
#                            demography = gg.demography,
#                            population = gg.population)
# 
# gg.habitat.dynamics <- habitat_dynamics()
# gg.demography.dynamics <- demography_dynamics(demo_dens_dep = demo_density_dependence(gg.trans.mat,
#                                                                                          fecundity_fraction = 0.8,
#                                                                                          survival_fraction = 0.8))
# gg.population.dynamics <- population_dynamics()
# gg.dynamics <- build_dynamics(habitat_dynamics = gg.habitat.dynamics,
#                                  demography_dynamics = gg.demography.dynamics,
#                                  population_dynamics = gg.population.dynamics,
#                                  order = c("habitat_dynamics",
#                                            "population_dynamics",
#                                            "demography_dynamics")
# )
# 
# sim_results <- simulation(state = gg.state,
#                           dynamics = gg.dynamics,
#                           timesteps = 100,
#                           replicates = 5)
# 
# ########################################################################

# ################# WITH DETERMINISTIC CHANGES TO HABITAT #####################
# 
# gg.habitat <- build_habitat(habitat_suitability = gg.hab.suit,
#                             carrying_capacity = gg.hab.k,
#                             misc = NULL)
# gg.demography <- build_demography(transition_matrix = gg.trans.mat,
#                                   scale = 'local', 
#                                   habitat_suitability = gg.hab.suit,
#                                   dispersal_parameters = gg.disp.param,
#                                   misc = NULL)
# gg.population <- build_population(population_raster = gg.pop)
# gg.state <- build_state(habitat = gg.habitat,
#                         demography = gg.demography,
#                         population = gg.population)
# 
# gg.habitat.dynamics <- habitat_dynamics()
# gg.demography.dynamics <- demography_dynamics()
# gg.population.dynamics <- population_dynamics()
# gg.dynamics <- build_dynamics(habitat_dynamics = gg.habitat.dynamics,
#                               demography_dynamics = gg.demography.dynamics,
#                               population_dynamics = gg.population.dynamics,
#                               order = c("habitat_dynamics",
#                                         "population_dynamics",
#                                         "demography_dynamics")
# )
# 
# sim_results <- simulation(state = gg.state,
#                           dynamics = gg.dynamics,
#                           timesteps = 100,
#                           replicates = 5)
# 
# ########################################################################
# 
# ################# WITH DETERMINISTIC CHANGES IN DEMOGRAPHICS #####################
# 
# gg.habitat <- build_habitat(habitat_suitability = gg.hab.suit,
#                             carrying_capacity = gg.hab.k,
#                             misc = NULL)
# gg.demography <- build_demography(transition_matrix = gg.trans.mat,
#                                   scale = 'local', 
#                                   habitat_suitability = gg.hab.suit,
#                                   dispersal_parameters = gg.disp.param,
#                                   misc = NULL)
# gg.population <- build_population(population_raster = gg.pop)
# gg.state <- build_state(habitat = gg.habitat,
#                         demography = gg.demography,
#                         population = gg.population)
# 
# gg.habitat.dynamics <- habitat_dynamics()
# gg.demography.dynamics <- demography_dynamics(demo_dens_dep = demo_density_dependence(gg.trans.mat,
#                                                                                       fecundity_fraction = 0.8,
#                                                                                       survival_fraction = 0.8))
# gg.population.dynamics <- population_dynamics()
# gg.dynamics <- build_dynamics(habitat_dynamics = gg.habitat.dynamics,
#                               demography_dynamics = gg.demography.dynamics,
#                               population_dynamics = gg.population.dynamics,
#                               order = c("habitat_dynamics",
#                                         "population_dynamics",
#                                         "demography_dynamics")
# )
# 
# sim_results <- simulation(state = gg.state,
#                           dynamics = gg.dynamics,
#                           timesteps = 100,
#                           replicates = 5)
# 
# ##################################################################################

plot(sim_results)

plot(sim_results, stage = 2)

plot(sim_results, stage = 0)

plot(sim_results[1], type = "raster", stage = 0, timesteps = c(10,20,30,40,50,60,70,80,90))

plot(sim_results[1], type = "raster", stage = 2)

plot(sim_results[1], type = "raster", stage = 3, animate = TRUE)

plot(sim_results[1], object = "habitat_suitability")

plot(sim_results[1], object = "carrying_capacity")

plot(sim_results[[1]][[10]]$population$population_raster)

