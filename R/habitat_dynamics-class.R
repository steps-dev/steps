#' Change the habitat in a state object
#'
#' @description A 'habitat dynamics' object is used to modify habitat suitability (or carrying capacity) of a landscape in space and time.
#' It is a sub-component of a \link[steps]{dynamics} object and is executed in each timestep of an experiment.
#'
#' @rdname habitat_dynamics
#'
#' @param habitat_dynamics_function A function that operates on a state object to change habitat at specified timesteps. User may enter a custom function or select a pre-defined module - see documentation. 
#' @param x an object to print or test as an habitat_dynamic object
#' @param ... further arguments passed to or from other methods
#' @param state a state object to apply the habitat function to
#' @param timestep the timestep in the experiment to apply the habitat function to the state object
#' @param habitat_suitability a raster layer or stack containing habitat suitability for each cell
#' @param disturbance_layers a raster stack with fire disturbances used to alter the habitat object in the experiment (number of layers must match the intended timesteps in the experiment)
#' @param effect_time the number of timesteps that the disturbance layer will act on the habitat object
#'
#' @return An object of class \code{habitat_dynamics}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' mat <- matrix(c(0.000,0.000,0.302,0.302,
#'                 0.940,0.000,0.000,0.000,
#'                 0.000,0.884,0.000,0.000,
#'                 0.000,0.000,0.793,0.793),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#' 
#' pop <- stack(replicate(4, ceiling(r * 0.2)))
#' 
#' dist <- r*0
#' 
#' dist[sampleRandom(dist, size=100, na.rm=TRUE, sp=TRUE)] <- 1
#' 
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"),
#'                               carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(transition_matrix = mat,
#'                                     dispersal_parameters = rlnorm(1))
#' test_population <- build_population(pop)
#' 
#' test_state <- build_state(test_habitat, test_demography, test_population)
#'
#' example_function <- function (state, timestep) {
#' state
#' }
#' 
#' no_habitat_dynamics <- as.habitat_dynamics(example_function)

as.habitat_dynamics <- function (habitat_dynamics_function) {
  stopifnot(inherits(habitat_dynamics_function,"function"))
  set_class(habitat_dynamics_function, "habitat_dynamics")
}

#' @rdname habitat_dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'habitat dynamics'
#'   
#' is.habitat_dynamics(no_habitat_dynamics)

is.habitat_dynamics <- function (x) {
  inherits(x, 'habitat_dynamics')
}

#' @rdname habitat_dynamics
#'
#' @export
#'
#' @examples
#'
#' print(no_habitat_dynamics)

print.habitat_dynamics <- function (x, ...) {
  cat("This is a habitat_dynamics object")
}

##########################
### internal functions ###
##########################




####################################
### pre-defined module functions ###
####################################

#' @rdname habitat_dynamics
#'
#' @export
#' 
#' @examples
#' 
#' # Use the no_habitat_dynamics object as a placeholder as it 
#' # does not modify the habitat object:
#' 
#' test_state2 <- no_habitat_dynamics(test_state, 1)
#' 
#' identical(test_state, test_state2)

no_habitat_dynamics <- function (state, timestep) {
  state
}

#' @rdname habitat_dynamics
#'
#' @export
#' 
#' @examples
#' 
#' # Use the fire_habitat_dynamics object to modify the  
#' # habitat using spatial fire history:
#' 
#' test_state2 <- fire_habitat_dynamics(habitat_suitability = r / cellStats(r, "max"),
#'                                     disturbance_layers = dist,
#'                                     effect_time = 1)

fire_habitat_dynamics <- function (habitat_suitability, disturbance_layers, effect_time=1) {
  
  habitat_dynamics <- function (state, timestep) {
    
    original_habitat <- habitat_suitability
    
    if (raster::nlayers(disturbance_layers) < timestep ) {
      stop("The number of disturbance layers must match the \nnumber of timesteps in the experiment")
    }

    modified_habitat <- original_habitat * raster::overlay(disturbance_layers[[utils::tail(seq_len(timestep), effect_time)]], fun=prod)
    names(modified_habitat) <- "Habitat"

    state$habitat$habitat_suitability <- modified_habitat
    
    state
    
  }
  
  as.habitat_dynamics(habitat_dynamics)
  
}
