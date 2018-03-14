#' Create a habitat object to use in a state object
#'
#' @description A habitat object is used to store spatially-explicit information on habitat suitability and the carrying_capacity of a landscape.
#' It is a sub-component of a \link[dhmpr]{state} object and is modified in each timestep of an experiment.
#' 
#' @rdname habitat
#' 
#' @param habitat_suitability A raster layer or stack containing habitat suitability for each cell
#' @param carrying_capacity A raster layer specifying carrying capacity values for each cell
#' @param x an object to print or test as an habitat object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{habitat}
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
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' # Modify raster to contain values between 0 and 1  
#' 
#' hs <- r / cellStats(r, "max")
#' 
#' # Modify raster to contain values for maximum population in each cell
#' 
#' k <- ceiling(r * 0.1)
#' 
#' # Construct the habitat object
#' 
#' test_habitat <- build_habitat(habitat_suitability = hs, carrying_capacity = k)

build_habitat <- function (habitat_suitability, carrying_capacity) {
  #INSERT CHECKS AND OBJECT TRANSFORMATIONS
  habitat <- list(habitat_suitability = habitat_suitability,
                  carrying_capacity = carrying_capacity)
  set_class(habitat, "habitat")
}

#' @rdname habitat
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'habitat'
#'   
#' is.habitat(test_habitat)

is.habitat <- function (x) {
  inherits(x, 'habitat')
}

#' @rdname habitat
#'
#' @export
#'
#' @examples
#' 
#' # Print information about the 'habitat' object
#'
#' print(test_habitat)

print.habitat <- function (x, ...) {
  cat("This is a habitat object")
}
