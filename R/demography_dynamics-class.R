#' Modify the demography in a state object.
#' 
#' A \code{demography_dynamics} object is used to modify life-stage transition
#' matrices - for example, adding stochasticity.
#' 
#' A \code{demography_dynamics} object is a sub-component of a \link[steps]{dynamics}
#' object and is executed in each timestep of a simulation. Note, some dynamics
#' functions can be executed at non-regular intervals (i.e. only timesteps
#' explicitly defined by the user). The \code{build_demography_dynamics} function is
#' used to construct a demography dynamics object consisting of several demographic
#' dynamics functions and their associated parameters. These functions specify how
#' the demography in the state object will be modified throughout a simulation.
#'
#' @rdname demography_dynamics
#'
#' @param ... Functions that operates on a state object to change demography
#' at specified timesteps. A user may enter custom functions or select
#' pre-defined modules - see examples. 
#' @param x A \code{demography_dynamics} object to print or test.
#'
#' @return An object of class \code{demography_dynamics}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' # Select existing dynamic functions to be run on the demography
#' # in a simulation and specify input parameters: 
#' env_stoch <- environmental_stochasticity(transition_matrix = egk_mat,
#'                                          stochasticity = 0.5)
#'                                               
#' dem_dens <- density_dependence(transition_matrix = egk_mat)
#' 
#' # Construct a demography dynamics object
#' dem_dynamics <- demography_dynamics(env_stoch, dem_dens)

demography_dynamics <- function (...) {
  
  dots <- list(...)
  
  # run checks on the functions passed in to make sure they are legit
  
  demo_dynamics <- function (state, timestep) {
    
    if (!is.null(unlist(dots))){
      for (fun in dots) {
        state <- fun(state, timestep)
      }
    }
    
    state
    
  }
  
  as.demography_dynamics(demo_dynamics)
  
}

as.demography_dynamics <- function (demography_dynamics_function) {
  as_class(demography_dynamics_function, "demography_dynamics", "function")
}

#' @rdname demography_dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'demography_dynamics'
#' is.demography_dynamics(test_demo_dynamics)

is.demography_dynamics <- function (x) {
  inherits(x, 'demography_dynamics')
}

#' @rdname demography_dynamics
#'
#' @export
#'
#' @examples
#' 
#' # Print details about the 'demography_dynamics' object
#' print(test_demo_dynamics)

print.demography_dynamics <- function (x, ...) {
  cat("This is a demography_dynamics object")
}