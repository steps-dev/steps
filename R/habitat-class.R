#' Create a habitat object to contain spatial information on habitat suitability and carrying capacity of a landscape
#'
#' @param habitat_suitability A raster layer or stack containing habitat suitability for each cell
#' @param carrying_capacity A raster layer specifying carrying capacity values for each cell
#'
#' @return An object of class \code{habitat}
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

build_habitat <- function (habitat_suitability, carrying_capacity) {
  habitat <- list(habitat_suitability = habitat_suitability,
                  carrying_capacity = carrying_capacity)
  set_class(habitat, "habitat")
}

#' Print details of a habitat object
#'
#' @param x an object to print or test as an habitat object
#' @param ... further arguments passed to or from other methods
#'
#' @export
#'
# @examples
# test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"), carrying_capacity = ceiling(r * 0.1))
# print(test_habitat)

print.habitat <- function (x, ...) {
  cat("This is a habitat object")
}