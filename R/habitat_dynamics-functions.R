#' Functions to modify the habitat in a landscape object.
#' 
#' Pre-defined functions to operate on habitat suitability (and carrying capacity if a function
#' is used) during a simulation.
#'
#' @name habitat_dynamics_functions
#' 
#' @seealso
#' \itemize{
#'   \item{\link[steps]{disturbance} to modify the suitability of a landscape with user provided
#'   spatially-explicit layers}
#'   \item{\link[steps]{fire_effects}}
#'   }
NULL

#' Disturbance
#' 
#' Modifies the landscape by multiplying habitat suitability values by a sum of previous
#' disturbances. Since disturbances can act in a single timestep, or have lasting effects,
#' the user can specify an 'effect time' of disturbances.
#'
#' @param disturbance_layers the name of spatial layer(s) in the landscape object with disturbances used
#'   to alter the habitat object for each timestep (number of layers must match the intended timesteps)
#' @param effect_time the number of timesteps that the disturbance layer will act on the habitat object
#'   (e.g. '3' will combine the effects of previous two timesteps to increase the overall effect) - the
#'   default is 1.
#'
#' @export
#' 
#' @examples
#'
#' # Road building (stored in the landscape object and called "roads") acts on the landscape
#' # each year.
#' 
#' \dontrun{
#' road_effect <- disturbance(disturbance_layers = "roads", effect_time = 1)
#' 
#' ls <- landscape(population = egk_pop, suitability = egk_hab, "roads" = egk_road)
#' 
#' pd <- population_dynamics(change = growth(egk_mat))
#' 
#' sim <- simulation(landscape = ls,
#'            population_dynamics = pd,
#'            habitat_dynamics = list(road_effect),
#'            timesteps = 20)
#'            
#' plot(sim, object = "suitability", type = "raster", timesteps = 1:9)
#' }

disturbance <- function (disturbance_layers, effect_time = 1) {
  
  dist_fun <- function (landscape, timestep) {
    
    if (raster::nlayers(landscape$suitability) > 1) {
      original_habitat <- landscape$suitability[[timestep]]
    } else {
      original_habitat <- landscape$suitability
    }
    
    if (raster::nlayers(landscape[[disturbance_layers]]) < timestep ) {
      stop("The number of disturbance layers must match the number of timesteps in the simulation")
    }
    
    # replace NA values with zeros
    landscape[[disturbance_layers]][is.na(landscape[[disturbance_layers]])] <- 0
    
    modified_habitat <- original_habitat * raster::overlay(landscape[[disturbance_layers]][[utils::tail(seq_len(timestep),
                                                                                                        effect_time)]],
                                                           fun = prod)
    names(modified_habitat) <- paste0("Habitat_", timestep)

    if (raster::nlayers(landscape$suitability) > 1) {
      landscape$suitability[[timestep]] <- modified_habitat
    } else {
      landscape$suitability <- modified_habitat
    }
    
    landscape
    
  }
  
  as.habitat_dynamics(dist_fun)
  
}

#' Fire effects with regeneration
#' 
#' Modifies the landscape by multiplying habitat suitability values by a weighted sum of previous
#' fire intensities based on a user specified regeneration function. By default, the regenerative
#' function is an inverse linear relationship to time, however, this function can be replaced with a response
#' that takes into account other factors of habitat restoration (e.g. growth/re-growth curves of vegetation).
#'
#' @param fire_layers the name(s) of spatial layer(s) in the landscape object with fire disturbances used
#'   to alter the habitat object for each timestep (number of layers must match the intended timesteps)
#' @param effect_time the number of timesteps that the fire layer will act on the habitat object
#' @param regeneration_function a function that determines how fast the landscape will regenerate after a
#'   fire event
#'
#' @export
#' 
#' @examples
#'
#' # Fire (stored in the landscape object and called "fires") acts on the landscape for
#' #five years with an exponentially decaying intensity.
#' 
#' \dontrun{
#' regen <- function (time) {-exp(time)}
#' 
#' plot(1:5, regen(1:5), type = "l")
#' 
#' fire <- fire_effects(fire_layers = "fires", effect_time = 5, regeneration_function = regen)
#' 
#' ls <- landscape(population = egk_pop, suitability = egk_hab, "fires" = egk_fire)
#' 
#' pd <- population_dynamics(change = growth(egk_mat))
#' 
#' sim <- simulation(landscape = ls,
#'            population_dynamics = pd,
#'            habitat_dynamics = list(fire),
#'            timesteps = 20)
#'            
#' plot(sim, object = "suitability", type = "raster", timesteps = 1:9)
#' }

fire_effects <- function (fire_layers,
                          effect_time = 3,
                          regeneration_function = function (time) {-time}) {
  
  dist_fun <- function (landscape, timestep) {
    
    if (raster::nlayers(landscape$suitability) > 1) {
      stop("This function only operates on landscape objects with a single initial habitat suitability layer -\n
           Please check that you have not provided a raster stack as a suitability component of a landscape object.")
    }
    
    if (is.null(steps_stash$orig_suitability)) {
      original_habitat <- steps_stash$orig_suitability <- landscape$suitability
    } else {
      original_habitat <- steps_stash$orig_suitability
    }
    
    if (raster::nlayers(landscape[[fire_layers]]) < timestep ) {
      stop("The number of fire layers must match the number of timesteps in the simulation")
    }
    
    # replace NA values with ones
    landscape[[fire_layers]][is.na(landscape[[fire_layers]])] <- 1
    
    # lags  
    years_since_fire <- c(0, seq_len(effect_time))
    
    # fire weights
    relative_regeneration <- regeneration_function(years_since_fire)
    regeneration <- rescale(relative_regeneration)
    
    lags <- timestep - years_since_fire
    valid <- lags > 0
    lags <- lags[valid]
    
    weights <- regeneration[valid]
    
    # apply the regeneration weights to previous fires
    fires_weighted <- (1 - landscape[[fire_layers]][[lags]]) * weights

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
  
  as.habitat_dynamics(dist_fun)
  
}


##########################
### internal functions ###
##########################

as.habitat_dynamics <- function (habitat_disturbance) {
  as_class(habitat_disturbance, "habitat_dynamics", "function")
}

rescale <- function (x) {
  x <- x - min(x)
  x / max(x)
}
