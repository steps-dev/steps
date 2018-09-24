#' Modify the habitat in a state object.
#'
#' A \code{habitat dynamics} object is used to modify habitat suitability
#' (or carrying capacity) of a landscape in space and time - for example,
#' by fire or intentional habitat modification.
#' 
#' The \code{build_habitat_dynamics} is a sub-component of a
#' \link[steps]{dynamics} object and is executed in each timestep of a
#' simulation. Note, some dynamics functions can be executed at
#' non-regular intervals (i.e. only timesteps explicitly defined by the user).
#' The \code{build_habitat_dynamics} function is used to construct a
#' habitat dynamics object consisting of several habitat
#' dynamics functions and their associated parameters. These functions
#' specify how the habitat in the state object will be modified
#' throughout a simulation.
#'
#' @rdname habitat_dynamics
#'
#' @param ... Functions that operate on a state object to change habitat
#' at specified timesteps. User may enter a custom function or select a
#' pre-defined module - see examples. 
#' @param x an \code{habitat_dynamics} object to print or test.
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
#' # Select existing habitat functions to be run on the demography
#' # in a simulation and specify input parameters: 
#' fires <- disturbance_fires(habitat_suitability = egk_hab,
#'                            disturbance_layers = egk_dist,
#'                            effect_time = 3)
#' 
#' # Construct a habitat dynamics object
#' hab_dynamics <- habitat_dynamics(fires)

habitat_dynamics <- function (...) {
  
  dots <- list(...)
  
  # run checks on the functions they've passed in, make sure they are legit
  
  hab_dynamics <- function (state, timestep) {
    
    if (!is.null(unlist(dots))){
      for (fun in dots) {
        state <- fun(state, timestep)
      }
    }
    
    state
    
  }
  
  as.habitat_dynamics(hab_dynamics)
  
}

as.habitat_dynamics <- function (habitat_dynamics_function) {
  as_class(habitat_dynamics_function, "habitat_dynamics", "function")
}

#' @rdname habitat_dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'habitat dynamics'
#' is.habitat_dynamics(test_hab_dynamics)

is.habitat_dynamics <- function (x) {
  inherits(x, 'habitat_dynamics')
}

#' @rdname habitat_dynamics
#'
#' @export
#'
#' @examples
#'
#' # Print details about the 'habitat_dynamics' object
#' print(test_hab_dynamics)

print.habitat_dynamics <- function (x, ...) {
  cat("This is a habitat_dynamics object")
}