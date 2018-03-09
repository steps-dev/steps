#' Change the demography in a state object
#'
#' @param demography_dynamics_function A function that operates on a state object to change demography at specified timesteps. User may enter a custom function or select a pre-defined module - see documentation. 
#'
#' @return An object of class \code{demography_dynamics}
#' @export
#'
#' @examples
#' 
#' library(raster)
#' library(dhmpr)
#'
#' example_function <- function (state, timestep) {
#' state
#' }
#' no_demography_dynamics <- as.demography_dynamics(example_function)

as.demography_dynamics <- function (demography_dynamics_function) {
  stopifnot(inherits(demography_dynamics_function,"function"))
  set_class(demography_dynamics_function, "demography_dynamics")
}

#' Print details of a demography_dynamics object
#'
#' @param x an object to print or test as an demography_dynamic object
#' @param ... further arguments passed to or from other methods
#'
#' @export
#'
# @examples
# example_function <- function (state, timestep) {
# state
# }
# no_demography_dynamics <- as.demography_dynamics(example_function)
# print(no_demography_dynamics)

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