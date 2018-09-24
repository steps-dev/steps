#' Modify the population in a state object.
#' 
#' A \code{population_dynamics} object is used to modify species populations
#' in space and time.
#' 
#' A population_dynamics object is a sub-component of a \link[steps]{dynamics}
#' object and is executed in each timestep of a simulation.  Note, some dynamics
#' functions can be executed at non-regular intervals (i.e. only timesteps
#' explicitly defined by the user). The \code{population_dynamics} function is
#' used to construct a population dynamics object consisting of several population
#' dynamics functions and their associated parameters. These functions specify how
#' the population in the state object will be modified throughout a simulation. 
#'
#' @rdname population_dynamics
#'
#' @param change a \link[steps]{population_dynamics_functions} to define how population growth occurs at each timestep
#' @param disp a function to define how the population disperses at each timestep (default is exponential kernel)
#' @param mod a function to define any deterministic changes to the population - such as translocation - at each timestep
#' @param dens_dep a function to control density dependence effects on the population at each timestep
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
#' uses default population growth function based on transition matrices
#' pop_dynamics <- population_dynamics()

population_dynamics <- function (change = simple_growth(),
                                 disp = NULL,
                                 mod = NULL,
                                 dens_dep = NULL) {
  
  pop_dynamics <- function (state, timestep) {
    
    state <- change(state, timestep)
    
    if (!is.null(disp))
      state <- disp(state, timestep)
    
    if (!is.null(mod))
      state <- mod(state, timestep)
    
    if (!is.null(dens_dep))
      state <- dens_dep(state, timestep)
    
    state
  }
  
  as.population_dynamics(pop_dynamics)
  
}

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
#' is.population_dynamics(test_pop_dynamics)

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
#' print(test_pop_dynamics)

print.population_dynamics <- function (x, ...) {
  cat("This is a population_dynamics object")
}