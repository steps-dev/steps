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
#' @param habitat A \link[steps]{habitat} object.
#' @param demography A \link[steps]{demography} object.
#' @param population A \link[steps]{population} object.
#' @param object A state object to print or test.
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
#' # Create a raster layer
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' # Create a life-stage matrix
#' mat <- matrix(c(0.000,0.000,0.302,0.302,
#'                 0.940,0.000,0.000,0.000,
#'                 0.000,0.884,0.000,0.000,
#'                 0.000,0.000,0.793,0.793),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#' 
#' # Create initial populations - count must match number of life-stages
#' pop <- stack(replicate(4, ceiling(r * 0.2)))
#' 
#' # Create habitat, demography, and population objects
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"),
#'                               carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(transition_matrix = mat)
#' test_population <- build_population(pop)
#' 
#' # Construct the state object
#' test_state <- build_state(test_habitat, test_demography, test_population)

build_state <- function (habitat, demography, population) {
  # CHECK OBJECT TYPES
  check_habitat_matches_population(habitat, population)
  check_demography_matches_population(demography, population)
  state <- list(habitat = habitat,
                demography = demography,
                population = population)
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

is.state <- function (object) {
  inherits(object, 'state')
}

#' @rdname state
#'
#' @export
#'
#' @examples
#' 
#' # Print information about the 'state' object
#' print(test_state)

print.state <- function (object, ...) {
  cat("This is a state object for a species with", ncol(object$demography$global_transition_matrix),
      "life stage(s) across a landscape of", lenght(object$habitat$habitat_suitability), "total cells.")
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
