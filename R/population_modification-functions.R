#' How the population is modified in a landscape.
#'
#' Pre-defined functions to define population modification (e.g. translocation) during a simulation.
#'
#' @name population_modification_functions
#'
#' @param stages which life-stages are modified - default is all
#' @param source_layer the name of a spatial layer in the landscape object with the locations
#'   and number of individuals to translocate from. Note, this layer will have only zero
#'   values if individuals are being introduced from outside the study area
#' @param sink_layer the name of a spatial layer in the landscape object with the locations
#'   and number of individuals to translocate to.
#' @param effect_timesteps which timesteps in a single simulation do the translocations
#'   take place
#'
#' @rdname population_modification_functions
#'
#' @export
#' 
#' @examples
#' 
#' # Use the translocation_population_dynamics object to modify the  
#' # population using translocations:
#' 
#' test_ca_dispersal <- translocation(source_layer = "egk_source",
#'                                        sink_layer = "egk_sink",
#'                                        stages = NULL,
#'                                        effect_timesteps = 1)

translocation <- function (source_layer, sink_layer, stages = NULL, effect_timesteps = 1) {
  
  pop_dynamics <- function (landscape, timestep) {
    
    if (timestep %in% effect_timesteps) {
      
      population_raster <- landscape$population
      nstages <- raster::nlayers(population_raster)
      
      # get population as a matrix
      idx <- which(!is.na(raster::getValues(population_raster[[1]])))
      population_matrix <- raster::extract(population_raster, idx)
      
      source <- raster::extract(landscape[[source_layer]], idx)
      sink <- raster::extract(landscape[[sink_layer]], idx)
      
      if (is.null(stages)) stages <- seq_len(nstages)
      
      for (stage in stages) {
        
        if (any(population_matrix[ , stage] - source < 0)) {
          warning("The proposed number of translocated individuals do not exist for\nlife-stage ", stage," in timestep ", timestep, " - only the maximum number of available\nindividuals in each cell will be translocated. Please check the\nspecified source and sink layers.")
        }
        
        population_matrix[ , stage] <- population_matrix[ , stage] + sink - pmin(source, population_matrix[ , stage])
        
      }
      
      # put back in the raster
      population_raster[idx] <- population_matrix
      
      landscape$population <- population_raster
      
    }
    
    landscape
    
  }
  
  as.population_translocation(pop_dynamics, info = list(source_layer = print(source_layer),
                                                        sink_layer = print(sink_layer),
                                                        stages = stages,
                                                        effect_timesteps = effect_timesteps))
  
}

##########################
### internal functions ###
##########################

as.population_translocation <- function (translocation, info = NULL) {
  as_class(translocation, "population_translocation", "function", info = info)
}

print.population_translocation <- function (x, ...) {
  print_info(x)
}