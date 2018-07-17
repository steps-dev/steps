#' Create a population object to use in a state object.
#' 
#' A population object contains information on how populations
#' (abundances of individuals) change in space and time.
#' 
#' A population object is a sub-component of a \link[steps]{state} object
#' and is modified in each timestep of a simulation. During a simulation,
#' a population object tracks changes in total individuals per grid cells
#' based on population dynamic functions selected or created by the user. 
#' 
#' @rdname population
#' 
#' @param population_raster A raster stack (grid cell-based) with one layer
#' for each life stage.
#' @param object A population object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{population}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' # Construct a raster object
#' r <- raster(system.file("external/test.grd",package="raster"))
#' 
#' # Create a stack of raster layers to represent each
#' # life-stage of a population structure (four in this case)
#' rs <- stack(replicate(4, r * 0.2))
#' 
#' # Construct the population object
#' test_population <- build_population(rs)

build_population <- function (population_raster) {
  population <- list(population_raster = population_raster)
  as.population(population)
}

#' @rdname population
#' 
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'population'
#' is.population(test_population)
 
is.population <- function (object) {
  inherits(object, 'population')
}

#' @rdname population
#'
#' @export
#'
#' @examples
#'
#' # Print information about the 'population' object
#' print(test_population)

print.population <- function (object, ...) {

  cat("This is a populaion object that contains ", nlayers(object[['population_raster']]),
      " life stage(s).")
  
}


##########################
### internal functions ###
##########################

as.population <- function (population) {
  as_class(population, "population", "list")
}