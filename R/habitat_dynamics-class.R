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

#' @export
no_habitat_dynamics <- function (state, timestep) {
  state
}

#' @export
fire_habitat_dynamics <- function (habitat_suitability, disturbance_layers, effect_time=1) {
  
  habitat_dynamics <- function (state, timestep) {
    
    original_habitat <- habitat_suitability
    
    if (nlayers(disturbance_layers) < timestep ) {
      stop("The number of disturbance layers must match the \nnumber of timesteps in the experiment")
    }

    modified_habitat <- original_habitat * overlay(disturbance_layers[[tail(seq_len(timestep), effect_time)]], fun=prod)
    names(modified_habitat) <- "Habitat"

    state$habitat$habitat_suitability <- modified_habitat
    
    state
    
  }
  
  as.habitat_dynamics(habitat_dynamics)
  
}
