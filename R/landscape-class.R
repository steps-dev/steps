#' Create a landscape object to use in a state object.
#'
#' A landscape object is used to store spatially-explicit information on landscape
#' suitability and the carrying_capacity of a landscape.
#' 
#' A landscape object is a sub-component of a \link[steps]{state} object and
#' is modified in each timestep of a simulation. During a simulation, a
#' landscape object tracks changes in landscape suitability or carrying capacity
#' based on landscape dynamic functions selected or created by the user.
#' 
#' @rdname landscape
#' 
#' @param population A raster stack (grid cell-based) with one layer
#' for each life stage.
#' @param suitability An optional raster layer or stack containing landscape
#'  suitability values for all cells in a landscape. Note, using a raster
#'  stack assumes that the user has provided a landscape layer for each
#'  intended timestep in a simulation.
#' @param carrying_capacity An optional raster layer specifying carrying capacity
#'  values for all cells in a landscape.
#' @param x A landscape object.
#' @param ... In the lanndscape function, named raster objects representing different
#'  aspects of the landscape used to modify the landscape object in a
#'  simulation. Note, this is intended to store objects that are accessed
#'  and used to modify the landscape with a custom developed landscape dynamic
#'  function. Also, further arguments passed to or from other methods.
#'
#' @return An object of class \code{landscape}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' # Construct the landscape object
#' hab <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)

landscape <- function (population, suitability = NULL, carrying_capacity = NULL, ...) {
  landscape <- list(population = population,
                    suitability = suitability,
                    carrying_capacity = carrying_capacity,
                    ...)
  as.landscape(landscape)
}

#' @rdname landscape
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'landscape'
#' is.landscape(test_habitat)

is.landscape <- function (x) {
  inherits(x, 'landscape')
}

#' @rdname landscape
#'
#' @export
#'
#' @examples
#' 
#' # Print information about the 'landscape' object
#' print(test_habitat)

print.landscape <- function (x, ...) {
  r.dims <- dim(x[['population']])
  r.res <- raster::res(x[['population']])
  
  cat("This is a landscape object that contains initial populations for", raster::nlayers(x[['population']]),
      "stage(s). Each grid cell is", r.res[1], "by", r.res[2], "map",
      "units (based on projection) in size and the grid cells are arranged in",
      r.dims[1], "rows and", r.dims[2], "columns.")
  
}

##########################
### internal functions ###
##########################

as.landscape <- function (landscape) {
  as_class(landscape, "landscape", "list")
}