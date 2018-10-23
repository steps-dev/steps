#' Functions to modify the habitat in a landscape object.
#' 
#' Pre-defined functions to operate on a habitat and carrying capacity
#' during a simulation.
#'
#' @name habitat_dynamics_functions
#'
#' @param disturbance_layers a raster stack with fire disturbances used to alter the habitat object in the experiment (number of layers must match the intended timesteps in the experiment)
#' @param effect_time the number of timesteps that the disturbance layer will act on the habitat object
#'
#' @examples
#' 
#' library(steps)
#' 
#' @rdname habitat_dynamics_functions
#'
#' @export
#' 
#' @examples
#' 
#' # Use the disturbance function to modify the habitat using spatial
#' # fire history layers (stored in the landscape object):
#' 
#' test <- disturbance(habitat_suitability = r / cellStats(r, "max"),
#'                                     disturbance_layers = "fires",
#'                                     effect_time = 1)

disturbance <- function (disturbance_layers, effect_time = 1) {
  
  dist_fun <- function (landscape, timestep) {
    
    if (timestep == 1) original_habitat <- landscape[["orig_suitability"]] <- landscape$suitability
    else original_habitat <- landscape[["orig_suitability"]]
    
    if (raster::nlayers(landscape[[disturbance_layers]]) < timestep ) {
      stop("The number of disturbance layers must match the number of timesteps in the experiment")
    }

    modified_habitat <- original_habitat * raster::overlay(landscape[[disturbance_layers]][[utils::tail(seq_len(timestep), effect_time)]], fun = prod)
    names(modified_habitat) <- "Habitat"

    landscape$suitability <- modified_habitat
    
    landscape
    
  }
  
  as.habitat_disturbance(dist_fun)
  
}

##########################
### internal functions ###
##########################

as.habitat_disturbance <- function (habitat_disturbance) {
  as_class(habitat_disturbance, "habitat_dynamics", "function")
}

# as.habitat_stochastic_disturbance <- function (habitat_stochastic_disturbance) {
#   as_class(habitat_stochastic_disturbance, "habitat_dynamics", "function")
# }