#' Create a proportion dispersing function
#'
#' A proportion dispersing function generates the proportions of species
#' that disperse from cells based on landscape features.
#' 
#' The default \code{set_proportion_dispersing} function and parameters returns
#' full dispersal for all life stages. Additional proportion dispersing functions
#' are provided in the software for the user to select, however, a user may also
#' provide a custom written proportion dispersing function. Please see the tutorial
#' vignette titled "Creating custom *steps* functions" for information on how to
#' write custom functions for use in simulations.
#' 
#' @name dispersal_proportion_function
#' @seealso
#' \itemize{
#'   \item{\link[steps]{set_proportion_dispersing} controls the proportions of each life-stage that disperse}
#'   \item{\link[steps]{density_dependence_dispersing} proportions of dispersing populations are controlled by
#'   approach to carrying capacity}
#'   }
NULL
#'  

#' Set proportions of populations dispersing
#'
#' This function allows a user to specify what proportions of populations in each life-stage disperse.
#' It operates similarly on all cells and in all timesteps throughout a simulation.
#'
#' @param proportions A single value or vector of proportions (between zero and one) of
#' individuals in each life stage that disperse - default is 1. If proportions are specified
#' as a single number, then all life-stages disperse with that proportion, however, a vector
#' of proportions (equal in length to the number of life-stages) can also be specified. To
#' prevent stages from dispersing, set corresponding values to zero.
#' 
#' @return An object of class \code{dispersal_proportion_function}
#' 
#' @export
#'
#' @examples
#' 
#' # Example of a proportion function that disperses no population in the first life stage,
#' # 50% of the second, and 90% of the 3rd. 
#' 
#' \dontrun{
#' prop_dispersal <- set_proportion_dispersing(proportions = c(0, 0.5, 0.9))
#' 
#' kb_dispersal <- kernel_dispersal(dispersal_proportion = prop_dispersal,
#'                       max_distance = 2000,
#'                       dispersal_kernel = exponential_dispersal_kernel(distance_decay = 1000))
#' 
#' ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)
#' 
#' pd <- population_dynamics(change = growth(egk_mat), dispersal = kb_dispersal)
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

set_proportion_dispersing <- function (proportions = 1) {
  
  disp_prop_fun <- function (landscape, timestep) {
   
    # get total life-stages
    n_stages <- raster::nlayers(landscape$population)

    dispersal_proportion <- int_or_proper_length_vector(proportions, n_stages, "proportions")

    dispersal_proportion
    
  }
  
  as.dispersal_proportion_function(disp_prop_fun)
}

#' Density-dependent proportions of populations dispersing
#' 
#' The proportion of populations dispersing will be density dependent in a simulation. Proportions
#' of populations in each life stage dispersing is adjusted based on available carrying capacity.
#' If life-stages are set by the \link[steps]{population_density_dependence_functions}, these
#' will be used to determine how close the population is to carrying capacity. If no
#' life-stages are set or density dependence is set to NULL in \link[steps]{population_dynamics},
#' the function will consider all life-stages in the calculation. 
#'
#' @param maximum_proportions A single value or vector of the maximum proportions (between zero
#' and one) of individuals in each life stage that disperse - default is 1. If maximum
#' proportions are specified as a single number, then all life-stages use that value, however,
#' a vector of maximum proportions (equal in length to the number of life-stages) can also be
#' specified. Maximum proportions are multiplied by the calculated proportions based on carrying
#' capacity so to prevent stages from dispersing, set corresponding values to zero.
#' 
#' @return An object of class \code{dispersal_proportion_function}
#' 
#' @export
#'
#' @examples
#' 
#' # Example of a proportion function that disperses all populations based on their approach
#' # to carrying capacity
#' 
#' \dontrun{
#' prop_dispersal <- density_dependence_dispersing()
#' 
#' kb_dispersal <- kernel_dispersal(dispersal_proportion = prop_dispersal,
#'                       max_distance = 2000,
#'                       dispersal_kernel = exponential_dispersal_kernel(distance_decay = 1000))
#' 
#' ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)
#' 
#' pd <- population_dynamics(change = growth(egk_mat), dispersal = kb_dispersal)
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

density_dependence_dispersing <- function (maximum_proportions = 1) {
  
  disp_prop_fun <- function (landscape, timestep) {

    # get total life-stages
    n_stages <- raster::nlayers(landscape$population)
    
    maximum_proportions <- int_or_proper_length_vector(maximum_proportions, n_stages, "maximum_proportions")
    
    # check for specified stages that contribute to density dependence
    if (!exists("density_dependence_stages")) {
      density_dependence_stages <- seq_len(n_stages)
    }
    
    # get non-NA cells
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))

    pop <- raster::getValues(landscape$population)
    
    # 22.01.20 - # cc <- get_carrying_capacity(landscape, timestep)
    # 22.01.20 - # cc <- raster::getValues(cc)
    cc <- raster::getValues(landscape$carrying_capacity) # 22.01.20

    dispersal_proportion <- rep(0, n_stages)
    
    for(i in density_dependence_stages) {
      proportion <- sum(pop[cell_idx, i]) / sum(cc[cell_idx])
      dispersal_proportion[i] <- tanh(proportion)
    }
    
    # check that maximum_proportions is a vector of the right length, with all values between 0 and 1
    dispersal_proportion * maximum_proportions
    
  }
  
  as.dispersal_proportion_function(disp_prop_fun)
}


##########################
### internal functions ###
##########################

as.dispersal_proportion_function <- function (dispersal_proportion_function) {
  as_class(dispersal_proportion_function, "dispersal_proportion_function", "function")
}
