#' Change the demography in a state object
#' 
#' @description A 'demography dynamics' object is used to modify life-stage transition matrices - adding stochasticity for example.
#' It is a sub-component of a \link[dhmpr]{dynamics} object and is executed in each timestep of an experiment.
#'
#' @rdname demography_dynamics
#'
#' @param demography_dynamics_function A function that operates on a state object to change demography at specified timesteps. User may enter a custom function or select a pre-defined module - see documentation. 
#' @param x an object to print or test as an demography_dynamic object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{demography_dynamics}
#' 
#' @export
#'
#' @examples
#' 
#' library(dhmpr)
#' library(raster)
#'
#' example_function <- function (state, timestep) {
#' state
#' }
#' 
#' no_demography_dynamics <- as.demography_dynamics(example_function)

as.demography_dynamics <- function (demography_dynamics_function) {
  stopifnot(inherits(demography_dynamics_function,"function"))
  set_class(demography_dynamics_function, "demography_dynamics")
}

#' @rdname demography_dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'demography dynamics'
#'   
#' is.demography_dynamics(no_demography_dynamics)

is.demography_dynamics <- function (x) {
  inherits(x, 'demography_dynamics')
}

#' @rdname demography_dynamics
#'
#' @export
#'
#' @examples
#' 
#' print(no_demography_dynamics)

print.demography_dynamics <- function (x, ...) {
  cat("This is a demography_dynamics object")
}

##########################
### internal functions ###
##########################




####################################
### pre-defined module functions ###
####################################

#' @export
no_demographic_dynamics <- function (state, timestep) {
  state
}

#' @export
envstoch_demographic_dynamics <- function (global_transition_matrix, stochasticity=0) {
  
  dim <- nrow(global_transition_matrix)
  idx <- which(global_transition_matrix != 0)
  recruitment_mask <- idx == ((dim ^ 2) - dim + 1)
  lower <- 0
  upper <- ifelse(recruitment_mask, 1, Inf)
  vals <- global_transition_matrix[idx]
  
  if (is.matrix(stochasticity)) {
    stopifnot(identical(dim(global_transition_matrix), dim(stochasticity)))
    stopifnot(identical(which(stochasticity != 0), idx))
    stochasticity <- stochasticity[idx]
  }
  
  demographic_dynamics <- function (state, timestep) {
    
    transition_matrix <- global_transition_matrix
    
    
    transition_matrix[idx] <- extraDistr::rtnorm(length(idx),
                                                 vals,
                                                 stochasticity,
                                                 a = lower,
                                                 b = upper)
    
    state$demography$transition_matrix <- transition_matrix
    
    state
    
  }
  
  as.demography_dynamics(demographic_dynamics)

}