#' How the population responds to density dependence in a landscape.
#'
#' Pre-defined or custom functions to define population density dependence (e.g. ceiling) during a simulation.
#'
#' @name population_density_dependence_functions
#'
#' @param stages which life-stages contribute to density dependence - default is all
#'
#' @rdname population_density_dependence_functions
#'
#' @export
#' 
#' @examples
#' 
#' test_pop_dd <- ceiling_density()

ceiling_density <- function (stages = NULL) {

  pop_dynamics <- function (landscape, timestep) {
    
    population_raster <- landscape$population

    # Get non-NA cells
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    
    cc <- get_carrying_capacity(landscape, timestep)

    if (is.null(cc)) {
      stop ("carrying capacity must be specified in the landscape object to use ceiling_density",
            call. = FALSE)
    }
    
    # get population as a matrix
    population_matrix <- raster::extract(population_raster, idx)
    carrying_capacity <- raster::extract(cc, idx)
    
    # get degree of overpopulation, and shrink accordingly
    if (!is.null(stages)) {
      overpopulation <- as.vector(carrying_capacity) / rowSums(cbind(population_matrix[ , stages], rep(0, length(idx))))
    } else {
      overpopulation <- as.vector(carrying_capacity) / rowSums(population_matrix)
    }
    
    overpopulation[is.nan(overpopulation)] <- 0
    overpopulation <- pmin(overpopulation, 1)
    population <- sweep(population_matrix, 1, overpopulation, "*")

    # get whole integers
    population <- round_pop(population)
    
    # put back in the raster
    population_raster[idx] <- population
    
    landscape$population <- population_raster
    
    landscape
  }

  result <- as.population_density_dependence(pop_dynamics)
  
  attr(result, "density_dependence_stages") <- stages
  return(result)
}


##########################
### internal functions ###
##########################

as.population_density_dependence <- function (density_dependence) {
  as_class(density_dependence, "population_density_dependence", "function")
}
