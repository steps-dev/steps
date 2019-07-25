#' Create a landscape object.
#'
#' A landscape object is used to store spatially-explicit information on population,
#' habitat suitability, carrying_capacity and miscellaneous landscape information.
#' 
#' A landscape object is modified in each timestep of a simulation. During a simulation,
#' population, habitat suitability or carrying capacity in a landscape object are changed 
#' based on dynamic functions selected or created by the user.
#' 
#' @rdname landscape
#' 
#' @param population a raster stack (grid cell-based) with one layer for each life
#'  stage.
#' @param suitability an optional raster layer or stack (multiple layers) containing
#'  habitat suitability values for all cells in a landscape. Note, using a raster
#'  stack assumes that the user has provided a layer for each intended timestep
#'  in a simulation.
#' @param carrying_capacity an optional raster layer specifying carrying capacity
#'  values for all cells in a landscape or a function defining how carrying capacity
#'  is determined by habitat suitability.
#' @param ... named raster objects representing different aspects of the landscape
#'  used to modify the landscape object in a simulation. Note, this is intended to
#'  store objects that are accessed by dynamic functions and used to modify the
#'  landscape in a simulation. Also, further arguments passed to or from other methods.
#'
#' @return An object of class \code{landscape}
#' 
#' @export
#'
#' @examples
#' 
#' # Example of setting up a landscape object.
#' 
#' \dontrun{
#' ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)
#' 
#' pd <- population_dynamics(change = growth(egk_mat))
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

landscape <- function (population, suitability = NULL, carrying_capacity = NULL, ...) {
  if(is.null(population)) stop("Initial population rasters must be provided for the landscape object.")
  
  if(!is.null(suitability)) {
    check_raster_matches_population(suitability, population)
    check_raster_na_matches(suitability, population)
  }
  
  if(!is.null(carrying_capacity) & identical(class(carrying_capacity)[1], "RasterLayer")) {
    check_raster_matches_population(carrying_capacity, population)
    check_raster_na_matches(carrying_capacity, population)
  }

  landscape <- list(population = population,
                    suitability = suitability,
                    carrying_capacity = carrying_capacity,
                    ...)
  
  as.landscape(landscape)
}

is.landscape <- function (x) {
  inherits(x, 'landscape')
}

# #' @rdname landscape
# #'
# #' @export
# #'
# #' @examples
# #'
# #' # Print information about the 'landscape' object
# #' print(ls)
# 
# print.landscape <- function (x, ...) {
# 
#   r.dims <- dim(x[['population']])
#   r.res <- raster::res(x[['population']])
# 
#   cat("This is a landscape object that contains initial populations for", raster::nlayers(x[['population']]),
#       "stage(s).\nEach grid cell is", r.res[1], "by", r.res[2], "map",
#       "units (based on projection) in size\nand a total of", prod(r.dims),
#       "grid cells are arranged in", r.dims[1], "rows and", r.dims[2], "columns.")
# 
# }

##########################
### internal functions ###
##########################

as.landscape <- function (landscape) {
  as_class(landscape, "landscape", "list")
}

check_raster_matches_population <- function (raster, population) {
  if (!identical(raster::res(raster), raster::res(population))) stop("Landscape rasters do not have matching resolution. ",
                                                                     "This must be corrected before running a simulation.")
  if (!identical(raster::extent(raster), raster::extent(population))) stop("Landscape rasters do not have matching extents. ",
                                                                           "This must be corrected before running a simulation.")
}

check_raster_na_matches <- function (raster, population) {
  
  pop_na <- which(is.na(raster::getValues(population[[1]])))
  suit_layers <- raster::nlayers(raster)
  
  if (suit_layers > 1) {
    for (i in seq_len(suit_layers)) {
      ras_na <- which(is.na(raster::getValues(raster[[i]])))
      if (!identical(ras_na, pop_na)) stop(paste0("Landscape suitability raster for timestep ",
                                           i,
                                           " does not have matching NA cells. ",
                                           "This must be corrected before running a simulation."))
    }
  } else { 
    ras_na <- which(is.na(raster::getValues(raster[[1]])))
    if (!identical(ras_na, pop_na)) stop("Landscape rasters do not have matching NA cells. ",
                                         "This must be corrected before running a simulation.")
  }
}