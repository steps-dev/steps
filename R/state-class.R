#' Create a state object
#'
#' A state object is a static representation of habitat, population,
#' and demography in a single timestep.
#' 
#' A state object is modified in each timestep of a simulation based on
#' specified functions in a \link[steps]{dynamics} object. When first
#' created, the state object will hold initial values for habitat, population,
#' and demography which may or may not change throughout a simulation.
#' 
#' @rdname state
#'
#' @param population A \link[steps]{population} object.
#' @param habitat A \link[steps]{habitat} object.
#' @param demography A \link[steps]{demography} object.
#' @param x A state object to print or test.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{state}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' # Create habitat, demography, and population objects
#' hab <- habitat(egk_hab, egk_k)
#' dem <- demography(egk_mat)
#' pop <- population(egk_pop)
#' 
#' # Construct the state object
#' state <- state(pop, hab, dem)

state <- function (population, habitat, demography) {
  # CHECK OBJECT TYPES
  check_habitat_matches_population(habitat, population)
  check_demography_matches_population(demography, population)
  state <- list(population = population,
                habitat = habitat,
                demography = demography)
  as.state(state)
}

#' @rdname state
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'state'
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
#' # Print information about the 'state' object
#' print(test_state)

print.state <- function (x, ...) {
  cat("This is a state object for a species with", ncol(x$demography$global_transition_matrix),
      "life stage(s) across a landscape of", length(x$habitat$habitat_suitability), "total cells.")
}

##########################
### internal functions ###
##########################

as.state <- function (state) {
  as_class(state, "state", "list")
}

check_habitat_matches_population <- function (habitat, population) {
  hab_ras <- habitat$habitat_suitability
  pop_ras <- population$population_raster
  stopifnot(identical(raster::res(hab_ras), raster::res(pop_ras)))
  stopifnot(identical(raster::extent(hab_ras), raster::extent(pop_ras)))
}

check_demography_matches_population <- function (demography, population) {
  stopifnot(ncol(demography$global_transition_matrix) ==
                      raster::nlayers(population$population_raster))
}
