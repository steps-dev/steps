#' Create a population object to use in a state object
#'
#' @param population_raster A raster stack with one layer for each life stage
#'
#' @return An object of class \code{population}
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
#' test_population <- build_population(stack(replicate(4, test_habitat$carrying_capacity * 0.2)))

build_population <- function (population_raster) {
  population <- list(population_raster = population_raster)
  set_class(population, "population")
}

#' Print details of a population object
#'
#' @param x an object to print or test as a population object
#' @param ... further arguments passed to or from other methods
#'
#' @export
#'
# @examples
# test_population <- build_population(stack(replicate(4, test_habitat$carrying_capacity * 0.2)))
# print(test_population)

print.population <- function (x, ...) {
  cat("This is a population object")
}

#' Verify a population object
#' @param x an object to print or test as a population object
#' 
#' @export
#' 
#' @examples
#' 
##' is.population(pops)
##' 
is.population <- function (x) {
  inherits(x, 'population')
}