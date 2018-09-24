
# Rename constructors to drop the 'build_',:
build_population_dynamics -> population_dynamics
habitat_dynamics -> habitat_dynamics
state -> landscape
pop_density_dependence -> population_cap


# Usage:

egk_landscape <- landscape(population = egk_pop, suitability = NULL, carrying_capacity = NULL)
# these are rasters, other named rasters can be passed via dots

egk_pop_dyn <- population_dynamics(growth(egk_mat),
                                   dispersal = NULL,
                                   modification = NULL,
                                   density_dependence = NULL)
# each of these components is a function, or a list of functions

# habitat dynamics is just a list of any functions that modify elements of the
# landscape (except the population)

test <- simulation(egk_landscape,
           population_dynamics = egk_pop_dyn,
           habitat_dynamics = list(),
           timesteps = 20,
           replicates = 3)

# no parallel flag, just use the future package

# examples:

# set up a landscape

# most simple simulation, just growth
simulation(landscape(egk_pop),
           population_dynamics(growth(egk_mat)))

# add density dependence (population_cap() requires the carrying_capacity raster
# to be provided )
test <- simulation(landscape(egk_pop,
                     carrying_capacity = egk_k),
           population_dynamics(growth(egk_mat),
                               density_dependence = population_cap()))


# add dispersal incorrectly: kernel_dispersal(arrival_probability = 'both')
# needs a suitabilty raster. Expect an error
simulation(landscape(egk_pop,
                     carrying_capacity = egk_k),
           population_dynamics(growth(egk_mat),
                               dispersal = kernel_dispersal(arrival_probability = "both"),
                               density_dependence = population_cap()))

# add dispersal
simulation(landscape(egk_pop,
                     suitability = egk_hab,
                     carrying_capacity = egk_k),
           population_dynamics(growth(egk_mat),
                               dispersal = kernel_dispersal(),
                               density_dependence = population_cap()))

# add the effects of fires on habitat suitability (since the population doesn't
# depend on habitat suitability, this has no effect on population, yet).
# disturbance_fires requires the habitat_suitability raster to be specified

simulation(landscape(egk_pop,
                     egk_hab,
                     egk_k),
           population_dynamics(growth(egk_mat),
                               dispersal = kernel_dispersal(),
                               density_dependence = population_cap()),
           habitat_dynamics = list(disturbance_fires(egk_hab, egk_dist)))
library(future)
plan(multiprocess)
test <- simulation(
  landscape(egk_pop,
            egk_hab,
            egk_k),
  population_dynamics(growth(egk_mat,
                             global_stochasticity = 0.05,
                             local_stochasticity = 0.01),
                      dispersal = kernel_dispersal(),
                      density_dependence = population_cap()),
  habitat_dynamics = list(disturbance_fires(egk_hab, egk_dist)),
  timesteps = 20,
  replicates = 8)

plot(test)



# have the leslie matrix depend on habitat suitability, and therefore be
# affected by the fires. suitability_affects() requires the habitat suitability
# raster
simulation(landscape(egk_pop,
                     egk_hab,
                     egk_k),
           population_dynamics(change = suitability_affects("fecundity"),
                               dispersal = kernel_dispersal(),
                               density_dependence = population_cap()
           ),
           habitat_dynamics = list(disturbance_fires(egk_hab, egk_dist)))

# have the leslie matrices depend on predefined rasters the transition matrix
# for each raster cell built once from rasters, the rasters named in
# vital_rate_rasters() must be named arguments additional to landscape() (passed
# via dots)
simulation(landscape(egk_pop,
                     fecundity = fec_raster,
                     survival = surv_raster),
           population_dynamics(
             change = vital_rate_rasters(fecundity = "fecundity",
                                         survival = "survival"),
             dispersal = kernel_dispersal(),
             density_dependence = population_cap()
           ))

# implicit population density dependence: fecundity parameters of the leslie
# matrix depend on habitat suitability, which depends on population size and
# carrying capacity. population_density_affects() requires the carrying_capacity
# raster, and modifies whatever raster is named (in this case
# habitat_suitability, but it could be a user-defined raster, like
# 'grass_cover')
simulation(landscape(egk_pop,
                     egk_hab,
                     egk_k),
           population_dynamics(
             change = suitability_affects("fecundity"),
             dispersal = kernel_dispersal(),
             density_dependence = population_cap()
           ),
           habitat_dynamics = list(
             population_density_affects("habitat_suitability")
           ))


# feed in survival & fecundity rasters

# build with survival/fecundity dependent on e.g. habitat suitability



population_dynamics(growth(egk_mat),
                    dispersal = kernel_dispersal(),
                    modification = NULL,
                    density_dependence = population_cap())

               
