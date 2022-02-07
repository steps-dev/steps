# # statics: 
# 
# habitat <- list(habitat_suitability = NULL,
#                 carrying_capacity = NULL)
# demography <- list(transition_matrix = NULL,
#                    dispersal_parameters = NULL)
# population <- list(population_raster = NULL)
# 
# state <- list(habitat = habitat,
#               demography = demography,
#               population = population)
# 
# # dynamics:
# 
# # changes the habitat (may depend on)
# habitat_dynamics <- function (state, timestep) {
#   state
# }
# demographic_dynamics <- function (state, timestep) {
#   state
# }
# population_dynamics <- function (state, timestep) {
#   state <- dispersal_function(state, timestep)
#   state <- population_change_function(state, timestep)
#   state
# }
# 
# dynamics <- list(habitat_dynamics = habitat_dynamics,
#                  demographic_dynamics = demographic_dynamics,
#                  population_dynamics = population_dynamics)



##### COPIED #######
experiment <- function (state, dynamics, timesteps = 100) {
  # check stuff
  timesteps <- seq_len(timesteps)
  output_states <- iterate_system(state, dynamics, timesteps)
  set_class(output_states, "experiment_results")
}

##### COPIED #######
iterate_system <- function (state, dynamics, timesteps) {
  
  output_states <- list()
  
  for (timestep in timesteps) {
    for (dynamic_function in dynamics) {
      state <- dynamic_function(state, timestep)
    }
    output_states[[timestep]] <- state
  }
  
  output_states
  
}

##### COPIED #######
set_class <- function (x, class) {
  class(x) <- c(class, class(x))
  x
}

##### COPIED #######
dispersal_matrix <- function (locations, distance_decay = 0.5) {
  D <- as.matrix(dist(locations))
  dispersal_matrix <- exp(-D / distance_decay)
  sums <- colSums(dispersal_matrix)
  dispersal_matrix <- sweep(dispersal_matrix, 2, sums, "/")
  dispersal_matrix
}

##### COPIED #######
fast_population_dynamics <- function (state, timestep) {
 
  population_raster <- state$population$population_raster
  dispersal_parameters <- state$demography$dispersal_parameters
  transition_matrix <- state$demography$transition_matrix

  # get population as a matrix
  idx <- which(!is.na(getValues(population_raster[[1]])))
  population <- extract(population_raster, idx)

  # do population change
  population <- population %*% transition_matrix
  
  # do dispersal
  locations <- raster::xyFromCell(population_raster, idx)
  resolution <- mean(res(population_raster))
  dispersal_decay <- dispersal_parameters * resolution
  
  dispersal <- dispersal_matrix(locations, dispersal_decay)
  population <- dispersal %*% population
  
  # put back in the raster
  population_raster[idx] <- population
  
  state$population$population_raster <- population_raster
  state
  
}
fast_population_dynamics <- set_class(fast_population_dynamics, "population_dynamics")

##### COPIED #######
# no changes to the habitat
no_habitat_dynamics <- function (state, timestep) {
  state
}
no_habitat_dynamics <- set_class(no_habitat_dynamics, "habitat_dynamics")

##### COPIED #######
# no changes to the demographics
no_demographic_dynamics <- function (state, timestep) {
  state
}
no_demographic_dynamics <- set_class(no_demographic_dynamics, "demographic_dynamics")

##### COPIED #######
check_habitat_matches_population <- function (habitat, population) {
  hab_ras <- habitat$habitat_suitability
  pop_ras <- population$population_raster
  stopifnot(identical(res(hab_ras), res(pop_ras)))
  stopifnot(identical(extent(hab_ras), extent(pop_ras)))
}

##### COPIED #######
check_demography_matches_population <- function (demography, population) {
  stopifnot(identical(ncol(demography$transition_matrix),
                      nlayers(population$population_raster)))
}

##### COPIED #######
build_state <- function (habitat, demography, population) {
  check_habitat_matches_population(habitat, population)
  check_demography_matches_population(demography, population)
  state <- list(habitat = habitat,
                demography = demography,
                population = population)
  set_class(state, "state")
}

##### COPIED #######
build_dynamics <- function (population_dynamics,
                            habitat_dynamics = no_habitat_dynamics,
                            demographic_dynamics = no_demographic_dynamics) {
  dynamics <- list(habitat_dynamics = habitat_dynamics,
                   demographic_dynamics = demographic_dynamics,
                   population_dynamics = population_dynamics)
  set_class(dynamics, "dynamics")
}

##### COPIED #######
build_habitat <- function (habitat_suitability, carrying_capacity) {
  habitat <- list(habitat_suitability = habitat_suitability,
                  carrying_capacity = carrying_capacity)
  set_class(habitat, "habitat")
}

##### COPIED #######
build_demography <- function (transition_matrix, dispersal_parameters) {
  demography <- list(transition_matrix = transition_matrix,
                     dispersal_parameters = dispersal_parameters)
  set_class(demography, "demography")
}

##### COPIED #######
build_population <- function (population_raster) {
  population <- list(population_raster = population_raster)
  set_class(population, "population")
}

##### COPIED #######
print.dynamics <- function(x, ...) {
  cat("dynamics object")
}

##### COPIED #######
print.population_dynamics <- function (x, ...) {
  cat("population_dynamics object")
}

##### COPIED #######
print.habitat_dynamics <- function (x, ...) {
  cat("habitat_dynamics object")
}

##### COPIED #######
print.demographic_dynamics <- function (x, ...) {
  cat("demographic_dynamics object")
}

##### COPIED #######
print.state <- function (x, ...) {
  cat("state object")
}

##### COPIED #######
print.habitat <- function (x, ...) {
  cat("habitat object")
}

##### COPIED #######
print.demography <- function (x, ...) {
  cat("demography object")
}

##### COPIED #######
print.population <- function (x, ...) {
  cat("population object")
}

##### COPIED #######
print.experiment_results <- function (x, ...) {
  cat("experiment results object, for", length(x), "timesteps")
}

# ~~~~~


##### COPIED #######
# fake a demographic transition matrix of the right dimensions
fake_transition_matrix <- function (n_stages) {
  
  survival <- runif(n_stages, 0.7, 0.9)
  growth <- runif(n_stages - 1, 0.5, 0.7)
  recruitment <- rlnorm(1)
  
  # base matrix 
  transition_matrix <- diag(n_stages)
  
  # add growth
  growth_idx <- which(row(transition_matrix) == col(transition_matrix) + 1)
  transition_matrix[growth_idx] <- growth
  
  # columns sum to 1
  sums <- colSums(transition_matrix)
  transition_matrix <- sweep(transition_matrix, 2, sums, FUN = "/")
  
  # apply survival
  transition_matrix <- sweep(transition_matrix, 2, survival, FUN = "*")
  
  # add recruitment
  transition_matrix[1, n_stages] <- recruitment
  
  transition_matrix
}

library (raster)
r <- raster(system.file("external/test.grd", package="raster"))

koala_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"),
                               carrying_capacity = ceiling(r * 0.1))
koala_demography <- build_demography(fake_transition_matrix(4),
                                     rlnorm(1))
koala_population <- build_population(stack(replicate(4, koala_habitat$carrying_capacity * 0.2)))
koala_state <- build_state(koala_habitat, koala_demography, koala_population)

fast_approximation <- build_dynamics(as.population_dynamics(fast_population_dynamics), 
                                     as.habitat_dynamics(no_habitat_dynamics), 
                                     as.demography_dynamics(no_demographic_dynamics)
                                     )

my_results <- experiment(koala_state,
                         fast_approximation,
                         timesteps = 10)

rasters <- lapply(my_results, function (state) state$population$population_raster[[1]])
plot(stack(rasters))
