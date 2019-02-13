#' Define population dynamics.
#' 
#' A \code{population_dynamics} object is used to describe how populations
#' change in space and time.
#' 
#' A population_dynamics object is passed to \link[steps]{simulation}
#' and defines how populations change between timesteps. Note, some dynamics
#' functions can be executed at non-regular intervals (i.e. only timesteps
#' explicitly defined by the user). The \code{population_dynamics} function is
#' used to construct a population dynamics object consisting of several population
#' dynamics functions and their associated parameters. These functions specify how
#' the population in the landscape object will be modified throughout a simulation. 
#'
#' @rdname population_dynamics
#'
#' @param change \link[steps]{population_change_functions} to define how population growth occurs at each timestep
#' @param dispersal a single or list of \link[steps]{population_dispersal_functions} to define how the population disperses at each timestep (default is exponential kernel)
#' @param modification a single or list of \link[steps]{population_modification_functions} to define any deterministic changes to the population - such as translocation - at each timestep
#' @param density_dependence a single or list of \link[steps]{population_density_dependence_functions} to control density dependence effects on the population at each timestep
#' @param x a population_dynamic object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{population_dynamics}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' # Construct a population dynamics object - note non-specified parameters
#' # uses default population growth function based on transition matrices
#' pop_dynamics <- population_dynamics()

population_dynamics <- function (change,
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL) {
  
  if (!is.null(density_dependence)) {
    
    # transfer information about the density dependence- contributing stages to the functions that need it
    dd_stages <- attr(density_dependence, "density_dependence_stages")
    
    if (!is.null(dispersal)) {
      # create an object called "density_dependence_stages" in the environment
      # in which dispersal() runs
      environment(dispersal)$density_dependence_stages <- dd_stages
    }
    
  }
  
  pop_dynamics <- function (landscape, timestep) {
    
    landscape <- change(landscape, timestep)
    
    if (!is.null(dispersal))
      landscape <- dispersal(landscape, timestep)
    
    if (!is.null(modification))
      landscape <- modification(landscape, timestep)
    
    if (!is.null(density_dependence))
      landscape <- density_dependence(landscape, timestep)
    
    landscape
  }
  
  as.population_dynamics(pop_dynamics)
  
}

################# Working ###################
# population_dynamics <- function (change,
#                                  dispersal = NULL,
#                                  modification = NULL,
#                                  density_dependence = NULL) {
# 
#   ordered_functions <- list()
#   
#   pop_dynamics <- function (landscape, timestep) {
#     
#     for (dynamic_function in ordered_functions) {
#       landscape <- dynamic_function(landscape, timestep)
#     }
#     
#     landscape <- change(landscape, timestep)
# 
#     landscape
#   }
# 
#   as.population_dynamics(pop_dynamics)
# 
# }
#########################################

as.population_dynamics <- function (population_dynamics_function) {
  as_class(population_dynamics_function, "population_dynamics", "function")
}


#' @rdname population_dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'population dynamics'
#' is.population_dynamics(pop_dynamics)

is.population_dynamics <- function (x) {
  inherits(x, 'population_dynamics')
}

#' @rdname population_dynamics
#'
#' @export
#'
#' @examples
#'
#' # Print details about the 'population_dynamics' object 
#' print(pop_dynamics)

print.population_dynamics <- function (x, ...) {
  cat("This is a population_dynamics object")
}