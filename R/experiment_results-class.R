#' Run an experiment to make spatially-explicit population projections
#'
#' @param state A state object - static habitat, population, and demography in a timestep
#' @param dynamics A dynamics object - modules that change habitat, population, and demography during and experiment
#' @param timesteps Number of timesteps used in the experiment
#' 
#' @return An object of class \code{experiment_results}
#' @export
#'
#' @examples
#' 
#' library(raster)
#' library(dhmpr)
#' 
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"), carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(fake_transition_matrix(4), rlnorm(1))
#' test_population <- build_population(stack(replicate(4, test_habitat$carrying_capacity * 0.2)))
#' test_state <- build_state(test_habitat, test_demography, test_population)
#' fast_approximation <- build_dynamics(no_habitat_dynamics, no_demographic_dynamics, fast_population_dynamics)
#' results <- experiment(test_state, fast_approximation, timesteps = 10)

experiment <- function (state, dynamics, timesteps = 100) {
  # check stuff
  timesteps <- seq_len(timesteps)
  output_states <- iterate_system(state, dynamics, timesteps)
  set_class(output_states, "experiment_results")
}

#' Print details of a experiment object
#'
#' @param x an object to print or test as an experiment_results object
#' @param ... further arguments passed to or from other methods
#'
#' @export
#'
# @examples
# results <- experiment(test_state, fast_approximation, timesteps = 10)
# print(results)

print.experiment_results <- function (x, ...) {
  cat("This is an experiment results object, for", length(x), "timesteps")
}

##########################
### internal functions ###
##########################

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