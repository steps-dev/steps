#' Functions to modify the habitat in a landscape object.
#' 
#' Pre-defined functions to operate on a habitat and carrying capacity
#' during a simulation.
#'
#' @name habitat_dynamics_functions
#'
#' @param disturbance_layers a raster stack with disturbances (e.g logging) used to alter the habitat object in the experiment (number of layers must match the intended timesteps in the experiment)
#' @param effect_time the number of timesteps that the disturbance layer will act on the habitat object
#' @param fire_layers a raster stack with fire disturbances used to alter the habitat object in the experiment (number of layers must match the intended timesteps in the experiment)
#' @param lag the number of timesteps that the fire layer will act on the habitat object
#' @param regeneration_function a function that determines how fast the landscape will regenerate after a fire event
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
#' # layers (stored in the landscape object):
#' 
#' logging <- disturbance(disturbance_layers = "logging",
#'                     effect_time = 1)

disturbance <- function (disturbance_layers, effect_time = 1) {
  
  dist_fun <- function (landscape, timestep) {
    
    if (timestep == 1) original_habitat <- steps_stash$orig_suitability <- landscape$suitability
    else original_habitat <- steps_stash$orig_suitability
    
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

#' @rdname habitat_dynamics_functions
#'
#' @export
#' 
#' @examples
#' 
#' # Use the fire_effects function to modify the habitat using spatial
#' # fire layers (stored in the landscape object) and a regeneration
#' # function:
#' 
#' fire <- fire_effects(fire_layers = "fires",
#'                     lag = 5,
#'                     regeneration_function = function (time) {-time})

fire_effects <- function (fire_layers,
                          lag = 3,
                          regeneration_function = function (time) {-time}) {
  
  dist_fun <- function (landscape, timestep) {
    
    if (timestep == 1) original_habitat <- steps_stash$orig_suitability <- landscape$suitability
    else original_habitat <- steps_stash$orig_suitability
    
    if (raster::nlayers(landscape[[fire_layers]]) < timestep ) {
      stop("The number of disturbance layers must match the number of timesteps in the experiment")
    }
    
    # lags  
    years_since_fire <- c(0, seq_len(lag))
    
    # fire weights
    relative_regeneration <- regeneration_function(years_since_fire)
    regeneration <- rescale(relative_regeneration)
    
    lags <- timestep - years_since_fire
    valid <- lags > 0
    lags <- lags[valid]
    
    weights <- regeneration[valid]
    
    # apply the regeneration weights to previous fires
    fires_weighted <- landscape[[fire_layers]][[lags]] * weights
    
    # get habitat reduction in all previous years, downweighted by
    # regeneration time
    annual_impact <- 1 - fires_weighted
    
    # get the cumulative impact
    if(raster::nlayers(annual_impact) == 1) total_impact <- annual_impact
    else total_impact <- prod(annual_impact)
    
    # apply the habitat reduction
    modified_habitat <- original_habitat * total_impact
    
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

rescale <- function (x) {
  x <- x - min(x)
  x / max(x)
}