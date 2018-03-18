#' Create a state object
#'
#' @description A state object represents a static representation of habitat, population, and demography in a timestep.
#' It is modified in each timestep of an experiment based on the specified \link[dhmpr]{dynamic} objects.
#' 
#' @rdname state
#'
#' @param habitat A habitat object - a habitat suitability raster layer or stack
#' @param demography A demography object - a stage-based transition matrix
#' @param population A population object - a raster stack with layers for each population stage
#' @param x an object to print or test as a state object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{state}
#' @export
#'
#' @examples
#' 
#' library(dhmpr)
#' library(raster)
#' 
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"), carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(transition_matrix = fake_transition_matrix(4), dispersal_parameters = rlnorm(1))
#' test_population <- build_population(stack(replicate(4, test_habitat$carrying_capacity * 0.2)))
#' 
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

#' @rdname state
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'state'
#'   
#' is.state(test_state)

is.state <- function (x) {
  inherits(x, 'state')
}

#' @rdname state
#'
#' @export
#'
#' @examples
#'
#' print(test_state)

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
  stopifnot(identical(ncol(demography$global_transition_matrix),
                      nlayers(population$population_raster)))
}