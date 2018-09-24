#' Create a dynamics object
#'
#' A dynamics object is a group of modular functions used modify habitat,
#' population, and demography in a timestep.
#' 
#' A dynamics object modifies a \link[steps]{state} object in each timestep
#' of a simulation based on specified habitat, demography, or population
#' dynamics functions.
#' 
#' @rdname dynamics
#' 
#' @param habitat_dynamics A module to alter the habitat object in a simulation
#' @param demography_dynamics A module to alter the habitat object in a simulation
#' @param population_dynamics A module to alter the population object in a simulation
#' @param order The order to apply the dynamics at each timestep in a simulation 
#' @param x an object to print or test as a dynamics object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{dynamics}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' dynamics <- dynamics(population_dynamics(),
#'                      habitat_dynamics(),
#'                      demography_dynamics())

dynamics <- function (population_dynamics,
                      habitat_dynamics,
                      demography_dynamics,
                      order = c("population_dynamics",
                                "habitat_dynamics",
                                "demography_dynamics")
) {
  #INSERT CHECKS FOR OBJECT CLASSES
  dynamics <- list(population_dynamics = population_dynamics,
                   habitat_dynamics = habitat_dynamics,
                   demography_dynamics = demography_dynamics)
  
  # get all the functions in a list, in the required order
  check_dynamics_order(order)
  dynamics <- lapply(order, get, envir = environment())
  as.dynamics(dynamics)
}

#' @rdname dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'dynamics'
#'   
#' is.dynamics(test_dynamics)

is.dynamics <- function (x) {
  inherits(x, 'dynamics')
}

#' @rdname dynamics
#'
#' @export
#'
#' @examples
#'
#' print(test_dynamics)

print.dynamics <- function(x, ...) {
  cat("This is a dynamics object")
}


##########################
### internal functions ###
##########################

as.dynamics <- function (dynamics) {
  as_class(dynamics, "dynamics", "list")
}

check_dynamics_order <- function (order) {
  sorted_order <- sort(order)
  expected <- c("demography_dynamics",
                "habitat_dynamics",
                "population_dynamics")
  if (!identical(sorted_order, expected)) {
    msg <- paste0("order must be a length-3 character vector giving the order ",
                  "in which to run the dynamic functions. It must contain each ",
                  "of the following strings once and only once:\n",
                  "'", paste(expected, collapse = "', '"), "'")
    stop (msg, call. = FALSE)
  }
}
