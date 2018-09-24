#' Create a habitat object to use in a state object.
#'
#' A habitat object is used to store spatially-explicit information on habitat
#' suitability and the carrying_capacity of a landscape.
#' 
#' A habitat object is a sub-component of a \link[steps]{state} object and
#' is modified in each timestep of a simulation. During a simulation, a
#' habitat object tracks changes in habitat suitability or carrying capacity
#' based on habitat dynamic functions selected or created by the user.
#' 
#' @rdname habitat
#' 
#' @param habitat_suitability A raster layer or stack containing habitat
#'  suitability values for all cells in a landscape. Note, using a raster
#'  stack assumes that the user has provided a habitat layer for each
#'  intended timestep in a simulation.
#' @param carrying_capacity A raster layer specifying carrying capacity
#'  values for all cells in a landscape.
#' @param misc Miscellaneous inputs used to modify the habitat object in a
#'  simulation. Note, this is intended to store objects that are accessed
#'  and used to modify the habitat with a custom developed habitat dynamic
#'  function.
#' @param x A habitat object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{habitat}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' # Construct the habitat object
#' hab <- habitat(habitat_suitability = egk_hab, carrying_capacity = egk_k)

habitat <- function (habitat_suitability, carrying_capacity, misc=NULL, ...) {
  habitat <- list(habitat_suitability = habitat_suitability,
                  carrying_capacity = carrying_capacity,
                  misc = misc)
  as.habitat(habitat)
}

#' @rdname habitat
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'habitat'
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
#' print(test_habitat)

print.habitat <- function (x, ...) {
  r.dims <- dim(x[['habitat_suitability']])
  r.res <- raster::res(x[['habitat_suitability']])
  
  cat("This is a habitat object that contains", raster::nlayers(x[['habitat_suitability']]),
      "habitat suitability layer(s). In both the habitat suitability and carrying",
      "capacity layers, each grid cell is", r.res[1], "by", r.res[2], "map",
      "units (based on projection) in size and the grid cells are arranged in",
      r.dims[1], "rows and", r.dims[2], "columns.")
  
}

##########################
### internal functions ###
##########################

as.habitat <- function (habitat) {
  as_class(habitat, "habitat", "list")
}