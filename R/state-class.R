#' Create a state object
#'
#' @description A state object represents a static representation of habitat, population, and demography in a timestep.
#' It is modified in each timestep of an experiment based on the specified dynamic objects.
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
#' library(steps)
#' library(raster)
#' 
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' mat <- matrix(c(0.000,0.000,0.302,0.302,
#'                 0.940,0.000,0.000,0.000,
#'                 0.000,0.884,0.000,0.000,
#'                 0.000,0.000,0.793,0.793),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#' 
#' pop <- stack(replicate(4, ceiling(r * 0.2)))
#' 
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"),
#'                               carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(transition_matrix = mat,
#'                                     dispersal_parameters = rlnorm(1))
#' test_population <- build_population(pop)
#' 
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
