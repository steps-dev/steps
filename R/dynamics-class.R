#' Create a dynamics object to run in an experiment
#'
#' @param population_dynamics A module to alter the population object in an experiment
#' @param habitat_dynamics A module to alter the habitat object in an experiment
#' @param demographic_dynamics A module to alter the habitat object in an experiment
#'
#' @return An object of class \code{dynamics}
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

build_dynamics <- function (habitat_dynamics,
                            demographic_dynamics,
                            population_dynamics
                            ) {
  #INSERT CHECKS FOR OBJECT CLASSES
  dynamics <- list(habitat_dynamics = habitat_dynamics,
                   demographic_dynamics = demographic_dynamics,
                   population_dynamics = population_dynamics)
  set_class(dynamics, "dynamics")
}

#' Print details of a dynamics object
#'
#' @param x an object to print or test as a dynamics object
#' @param ... further arguments passed to or from other methods
#'
#' @export
#'
# @examples
# fast_approximation <- build_dynamics(no_habitat_dynamics, no_demographic_dynamics, fast_population_dynamics)
# print(fast_approximation)

print.dynamics <- function(x, ...) {
  cat("This is a dynamics object")
}