#' How the population responds to density dependence in a landscape.
#'
#' Pre-defined functions to define population density dependence (e.g. cap) during a simulation.
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
#' test_pop_dd <- population_cap()

population_cap <- function (stages = NULL) {
  
  pop_dynamics <- function (landscape, timestep) {
    
    population_raster <- landscape$population

    # Get non-NA cells
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    
    if (exists("carrying_capacity_function", envir = steps_stash)) {
      landscape$carrying_capacity[idx] <- steps_stash$carrying_capacity_function(landscape$suitability[idx])
    }

    # get population as a matrix
    population_matrix <- raster::extract(population_raster, idx)
    carrying_capacity <- raster::extract(landscape$carrying_capacity, idx)
    
    # if (is.null(carrying_capacity)) {
    #   stop ("carrying capacity must be specified",
    #         call. = FALSE)
    # }
    
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
  
  as.population_density_dependence(pop_dynamics)
  
}

# IN PROGRESS...
# @rdname population_density_dependence_functions
#
# @export
# 
# @examples
# 
# test_pop_dd <- contest()

# contest <- function (transition_matrix,
#                      max_growth_rate = 1,
#                      stages = NULL) {
#   
#   pop_dynamics <- function (landscape, timestep) {
#     
#     population_raster <- landscape$population
#     
#     # Get non-NA cells
#     idx <- which(!is.na(raster::getValues(population_raster[[1]])))
#     
#     if (exists("carrying_capacity_function", envir = steps_stash)) {
#       landscape$carrying_capacity[idx] <- steps_stash$carrying_capacity_function(landscape$suitability[idx])
#     }
#     
#     # get population as a matrix
#     population_matrix <- raster::extract(population_raster, idx)
#     carrying_capacity <- raster::extract(landscape$carrying_capacity, idx)
# 
#     # get growth rates
#     if (!is.null(stages)) {
#       rate <- (max_growth_rate * as.vector(carrying_capacity)) / 
#         (max_growth_rate * rowSums(cbind(population_matrix[ , stages], rep(0, length(idx)))) - rowSums(cbind(population_matrix[ , stages], rep(0, length(idx)))) + as.vector(carrying_capacity)
#     } else {
#       rate <- max_growth_rate * as.vector(carrying_capacity) / max_growth_rate * rowSums(population_matrix) - rowSums(population_matrix) + as.vector(carrying_capacity)
#     }
#     
#     # get whole integers
#     population <- round_pop(population)
#     
#     # put back in the raster
#     population_raster[idx] <- population
#     
#     landscape$population <- population_raster
#     
#     landscape
#   }
# 
# as.population_density_dependence(pop_dynamics)
# 
# }

##########################
### internal functions ###
##########################

as.population_density_dependence <- function (density_dependence) {
  as_class(density_dependence, "population_dynamics", "function")
}
