#' Create a landscape object.
#'
#' A landscape object is used to store spatially-explicit information on population,
#' habitat suitability and carrying_capacity.
#' 
#' A landscape object is modified in each timestep of a simulation. During a simulation, a
#' landscape object tracks changes in population, habitat suitability or carrying capacity
#' based on dynamic functions selected or created by the user.
#' 
#' @rdname landscape
#' 
#' @param population A raster stack (grid cell-based) with one layer
#'  for each life stage.
#' @param suitability An optional raster layer or stack containing habitat
#'  suitability values for all cells in a landscape. Note, using a raster
#'  stack assumes that the user has provided a layer for each intended timestep
#'  in a simulation.
#' @param carrying_capacity An optional raster layer specifying carrying capacity
#'  values for all cells in a landscape or a function defining how carrying capacity
#'  is determined by habitat suitability.
#' @param x A landscape object.
#' @param ... Named raster objects representing different
#'  aspects of the landscape used to modify the landscape object in a
#'  simulation. Note, this is intended to store objects that are accessed
#'  and used to modify the landscape with custom developed dynamic
#'  functions. Also, further arguments passed to or from other methods.
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
#' egk_ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)

landscape <- function (population, suitability = NULL, carrying_capacity = NULL, ...) {
  if(is.null(population)) stop("Initial population rasters must be provided for the landscape object.")
  if(!is.null(suitability)) {
    check_raster_matches_population(suitability, population)
  }
  if(!is.null(carrying_capacity) & identical(class(carrying_capacity)[1], "RasterLayer")) {
    check_raster_matches_population(carrying_capacity, population)
  }
  if(!is.null(carrying_capacity) & identical(class(carrying_capacity)[1], "function")){
    assign("carrying_capacity_function", carrying_capacity, steps_stash)
    if(is.null(suitability)) stop("A carrying capacity function requires a suitability layer in the landscape object.")
    cell_idx <- which(!is.na(raster::getValues(suitability)))
    carrying_capacity <- suitability
    names(carrying_capacity) <- "k"
    carrying_capacity[cell_idx] <- steps_stash$carrying_capacity_function(suitability[cell_idx])
  }
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
#' is.landscape(egk_ls)

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
#' print(egk_ls)

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

check_raster_matches_population <- function (raster, population) {
  stopifnot(identical(raster::res(raster), raster::res(population)))
  stopifnot(identical(raster::extent(raster), raster::extent(population)))
}