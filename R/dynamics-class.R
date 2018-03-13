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
                            demography_dynamics,
                            population_dynamics,
                            order = c("demography_dynamics",
                                      "habitat_dynamics",
                                      "population_dynamics")
                            ) {
  #INSERT CHECKS FOR OBJECT CLASSES
  dynamics <- list(habitat_dynamics = habitat_dynamics,
                   demography_dynamics = demography_dynamics,
                   population_dynamics = population_dynamics)

  # get all the functions in a list, in the required order
  check_dynamics_order(order)
  dynamics <- lapply(order, get, envir = environment())
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


##########################
### internal functions ###
##########################

check_dynamics_order <- function (order) {
  sorted_order <- sort(order)
  expected <- c("demography_dynamics",
                "habitat_dynamics",
                "population_dynamics")
  if (!identical(sorted_order, expected)) {
    msg <- paste0("order must be a length-4 character vector giving the order ",
                  "in which to run the dynamic functions. It must contain each ",
                  "of the following strings once and only once:\n",
                  "'", paste(expected, collapse = "', '"), "'")
    stop (msg, call. = FALSE)
  }
}
