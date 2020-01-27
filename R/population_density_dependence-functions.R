#' How the population responds to density dependence in a landscape.
#'
#' Pre-defined or custom functions to define population density dependence (e.g. ceiling)
#' during a simulation. Please see the tutorial vignette titled "Creating custom *steps*
#' functions" for information on how to write custom functions for use in simulations.
#'
#' @name population_density_dependence_functions
#' 
#' @seealso
#' \itemize{
#'   \item{\link[steps]{ceiling_density} to cap populations at carrying capacities}
#'   }
NULL

#' Ceiling-based density dependence
#' 
#' In-built density dependence function that constrains the number of individuals in a cell
#' based on the carrying capacity of that cell in a timestep. Note, carrying_capacity must
#' be provided in the landscape object to use this function (see \link[steps]{landscape}).
#' Only specified stages that contribute to density dependence are considered in the
#' calculations and excess individuals are removed from only the contributing stages. This
#' type of density dependence only affects the population once it reaches the carrying
#' capacity. While population size is below carrying capacity, the population grows according
#' to the transition matrix. 
#'
#' @param stages which life-stages contribute to density dependence and are removed in a timestep
#'  - default is all
#'
#' @export
#' 
#' @examples
#' 
#' # Cap the population at carrying capacity with only the second and third
#' # life stage used in calculations to determine density dependence. 
#' 
#' \dontrun{
#' cap_population <- ceiling_density(stages = c(2, 3))
#' 
#' ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)
#' 
#' pd <- population_dynamics(change = growth(egk_mat), density_dependence = cap_population)
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

ceiling_density <- function (stages = NULL) {

  pop_dynamics <- function (landscape, timestep) {
 
    #browser()
    
    population_raster <- landscape$population

    # Get non-NA cells
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    
    # 22.01.20 - # cc <- get_carrying_capacity(landscape, timestep)
    cc <- landscape$carrying_capacity # 22.01.20

    if (is.null(cc)) {
      stop ("carrying capacity must be specified in the landscape object to use ceiling_density",
            call. = FALSE)
    }

    # get population as a matrix
    population_matrix <- raster::extract(population_raster, idx)
    carrying_capacity <- raster::extract(cc, idx)
    
    # get degree of overpopulation, and shrink accordingly
    if (is.null(stages))
      stages <- seq_len(ncol(population_matrix))
    
    overpopulation <- as.vector(carrying_capacity) / rowSums(population_matrix[ , stages, drop = FALSE])
    
    overpopulation[is.nan(overpopulation)] <- 0
    overpopulation <- pmin(overpopulation, 1)
    population_matrix[, stages] <- sweep(population_matrix[, stages, drop = FALSE], 1, overpopulation, "*")

    # get whole integers
    population_matrix <- round_pop(population_matrix)
    
    # put back in the raster
    population_raster[idx] <- population_matrix
    
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
