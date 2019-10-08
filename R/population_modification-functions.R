#' How the population is modified in a landscape.
#'
#' Pre-defined functions to define population modification (e.g. translocation) during a simulation.
#'
#' @name population_modification_functions
#' 
#' @seealso
#' \itemize{
#'   \item{\link[steps]{translocation} for specifying explicit spatial and temporal movements
#'   of populations}
#'   \item{\link[steps]{mortality} for specifying explicit spatial and temporal changes to populations}
#'   }
NULL

#' Translocate populations
#'
#' This function is used to move or introduce populations throughout a simulation. A user can
#' specify which life-stages will be affected (default is all) and in which timesteps the
#' translocations will take place. A warning will be generated if populations are not available
#' where specified to translocate from.
#'
#' @param origins_layer the name of a spatial layer in the landscape object with the locations
#'   and number of individuals to translocate from. Note, this layer will have only zero
#'   values if individuals are being introduced from outside the study area
#' @param destinations_layer the name of a spatial layer in the landscape object with the locations
#'   and number of individuals to translocate to. Note, this layer will have only zero
#'   values if individuals are being controlled (e.g. culling)
#' @param stages which life-stages are modified - default is all
#' @param effect_timesteps which timesteps in a single simulation do the translocations
#'   take place

#'
#' @export
#' 
#' @examples
#' 
#' # Modify populations in all life-stages using explicit layers of origin and destination populations
#' # in timesteps 5, 10, and 15.
#' 
#' \dontrun{
#' trans_pop <- translocation(origins_layer = "origins",
#'                            destinations_layer = "destinations",
#'                            stages = NULL,
#'                            effect_timesteps = c(5, 10, 15))
#' 
#' ls <- landscape(population = egk_pop,
#'                 suitability = NULL,
#'                 carrying_capacity = NULL,
#'                 "origins" = egk_origins,
#'                 "destinations" = egk_destinations)
#' 
#' pd <- population_dynamics(change = growth(egk_mat), modification = trans_pop)
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

translocation <- function (origins_layer, destinations_layer, stages = NULL, effect_timesteps = 1) {
  
  pop_dynamics <- function (landscape, timestep) {
    
    if (timestep %in% effect_timesteps) {
      
      population_raster <- landscape$population
      nstages <- raster::nlayers(population_raster)
      
      # get population as a matrix
      idx <- which(!is.na(raster::getValues(population_raster[[1]])))
      population_matrix <- raster::extract(population_raster, idx)
      
      origins <- raster::extract(landscape[[origins_layer]], idx)
      destinations <- raster::extract(landscape[[destinations_layer]], idx)
      
      if (is.null(stages)) stages <- seq_len(nstages)
      
      for (stage in stages) {
        
        warn_once(any(population_matrix[ , stage] - origins < 0),
                  paste("The proposed number of translocated individuals do not exist for\nlife-stage",
                        stage,
                        "- only the maximum number of available\nindividuals in each cell will be translocated. Please check the\nspecified origins and destination layers."),
                  warning_name = paste0("translocated_individuals_", stage))

        population_matrix[ , stage] <- population_matrix[ , stage] + destinations - pmin(origins, population_matrix[ , stage])
        
      }
      
      # put back in the raster
      population_raster[idx] <- population_matrix
      
      landscape$population <- population_raster
      
    }
    
    landscape
    
  }
  
  as.population_modification(pop_dynamics)
  
}


#' Directly affect populations
#' 
#' This function modifies a population by a mortality spatial layer included in a
#' steps landscape object. The mortality layer consists of values from 0???1 and
#' modifies the population by multiplying the population of a cell by the value of
#' the corresponding cell in a mortality layer. For example, a cell with ten
#' individuals before the mortality function is applied, and corresponding mortality
#' layer cell with a value of 0.2, would have two individuals remaining after
#' modification. Note, rounding also occurs after modification using a ceiling method
#' (i.e the largest whole integer is retained).

#' @param mortality_layer the name of spatial layer(s) in the landscape object with
#'   mortality proportions used to alter the populations for each timestep (number of
#'   layers must match the intended timesteps)
#' @param stages which life-stages are modified - default is all
#'
#' @export
#' 
#' @examples
#' # Modify populations in all life-stages with fire intensity.
#' 
#' \dontrun{
#' fire_mortal <- mortality(mortality_layer = "fire", stages = NULL)
#' 
#' ls <- landscape(population = egk_pop,
#'                 suitability = egk_hab,
#'                 carrying_capacity = egk_k,
#'                 "fire" = egk_fire)
#' 
#' pd <- population_dynamics(change = growth(egk_mat), modification = fire_mortal)
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

mortality <- function (mortality_layer, stages = NULL) {
  
  pop_dynamics <- function (landscape, timestep) {

      population_raster <- landscape$population
      nstages <- raster::nlayers(population_raster)
      
      # get population as a matrix
      idx <- which(!is.na(raster::getValues(population_raster[[1]])))
      population_matrix <- raster::extract(population_raster, idx)
      
      mortality_prop <- raster::extract(landscape[[mortality_layer]][[timestep]], idx)
      
      if (is.null(stages)) stages <- seq_len(nstages)
      
      for (stage in stages) {
        
        population_matrix[ , stage] <- ceiling(population_matrix[ , stage] * mortality_prop)
        
      }
      
      # put back in the raster
      population_raster[idx] <- population_matrix
      
      landscape$population <- population_raster

    
    landscape
    
  }
  
  as.population_modification(pop_dynamics)
  
}

##########################
### internal functions ###
##########################

as.population_modification <- function (modification) {
  as_class(modification, "population_modification", "function")
}