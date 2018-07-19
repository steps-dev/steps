#' Functions to modify the habitat in a state object.
#' 
#' Pre-defined functions to operate on a habitat and carrying capacity
#' during a simulation.
#'
#' @name habitat_dynamics_functions
#'
#' @param habitat_suitability a raster layer or stack containing habitat suitability for each cell
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
#' # Use the disturbance_fires function to modify the habitat using spatial
#' # fire history layers:
#' 
#' test_fires <- disturbance_fires(habitat_suitability = r / cellStats(r, "max"),
#'                                     disturbance_layers = dist,
#'                                     effect_time = 1)

disturbance_fires <- function (habitat_suitability, disturbance_layers, effect_time=1) {
  
  dist_fire_fun <- function (state, timestep) {
    
    original_habitat <- habitat_suitability
    
    if (raster::nlayers(disturbance_layers) < timestep ) {
      stop("The number of disturbance layers must match the \nnumber of timesteps in the experiment")
    }

    modified_habitat <- original_habitat * raster::overlay(disturbance_layers[[utils::tail(seq_len(timestep), effect_time)]], fun=prod)
    names(modified_habitat) <- "Habitat"

    state$habitat$habitat_suitability <- modified_habitat
    
    state
    
  }
  
  as.habitat_disturbance(dist_fire_fun)
  
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