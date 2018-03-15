#' Create a population object to use in a state object
#' 
#' @description A population object is used to store spatially-explicit information on species populations.
#' It is a sub-component of a \link[dhmpr]{state} object and is modified in each timestep of an experiment.
#' 
#' @rdname population
#' 
#' @param population_raster A raster stack with one layer for each life stage
#' @param x a population object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{population}
#' 
#' @export
#'
#' @examples
#' 
#' library(dhmpr)
#' library(raster)
#' 
#' # Construct a raster object
#' 
#' r <- raster(system.file("external/test.grd",package="raster"))
#' 
#' # Create a stack of raster layers to represent each
#' # life-stage of a population structure (four in this case)
#' 
#' rs <- stack(replicate(4, r * 0.2))
#' 
#' # Construct the population object
#' 
#' test_population <- build_population(rs)

build_population <- function (population_raster) {
  #ADD CHECKS AND OBJECT CONVERSIONS
  population <- list(population_raster = population_raster)
  set_class(population, "population")
}

#' @rdname population
#' 
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'population'
#' 
#' is.population(test_population)
 
is.population <- function (x) {
  inherits(x, 'population')
}

#' @rdname population
#'
#' @export
#'
#' @examples
#'
#' # Print information about the 'population' object
#' 
#' print(test_population)

print.population <- function (x, ...) {
  cat("This is a population object")
}