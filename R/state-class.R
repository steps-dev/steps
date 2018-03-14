#' Create a state object to represent a static representation of habitat, population, and demography in a timestep
#'
#' @param habitat A habitat object - a habitat suitability raster layer or stack
#' @param demography A demography object - a stage-based transition matrix
#' @param population A population object - a raster stack with layers for each population stage
#'
#' @return An object of class \code{state}
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

build_state <- function (habitat, demography, population) {
  # CHECK OBJECT TYPES
  check_habitat_matches_population(habitat, population)
  check_demography_matches_population(demography, population)
  state <- list(habitat = habitat,
                demography = demography,
                population = population)
  set_class(state, "state")
}

#' Print details of a state object
#'
#' @param x an object to print or test as a state object
#' @param ... further arguments passed to or from other methods
#'
#' @export
#'
# @examples
# test_state <- build_state(test_habitat, test_demography, test_population)
# print(test_state)

print.state <- function (x, ...) {
  cat("This is a state object")
}

##########################
### internal functions ###
##########################

check_habitat_matches_population <- function (habitat, population) {
  hab_ras <- habitat$habitat_suitability
  pop_ras <- population$population_raster
  stopifnot(identical(res(hab_ras), res(pop_ras)))
  stopifnot(identical(extent(hab_ras), extent(pop_ras)))
}

check_demography_matches_population <- function (demography, population) {
  stopifnot(identical(ncol(demography$transition_matrix),
                      nlayers(population$population_raster)))
}