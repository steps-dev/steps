#' Define population dynamics.
#' 
#' A \code{population_dynamics} object is used to describe how populations
#' change in space and time.
#' 
#' A population_dynamics object is passed to \link[steps]{simulation}
#' and defines how populations change between timesteps. Note, some dynamics
#' functions can be executed at non-regular intervals (i.e. only timesteps
#' explicitly defined by the user). The \code{population_dynamics} function is
#' used to construct an object with several population dynamics functions and
#' their associated parameters. These functions specify how the population in
#' the landscape object will be modified throughout a simulation. The dynamics
#' can be executed in any order that is specified by the user. It is cautioned
#' that the order of dynamics will have implications depending on whether the
#' user has assumed a post-breeding or pre-breeding census in the transition
#' matrix. For more information on this, please refer to Kendall et al, (2019)
#' \emph{Ecological Applications}.
#'
#' @rdname population_dynamics
#'
#' @param change \link[steps]{population_change_functions} to define how population
#'  growth occurs at each timestep
#' @param dispersal \link[steps]{population_dispersal_functions} to define how the
#'  population disperses at each timestep
#' @param modification \link[steps]{population_modification_functions} to define any
#'  deterministic changes to the population - such as translocations or population
#'  control - at each timestep
#' @param density_dependence \link[steps]{population_density_dependence_functions}
#'  to control density dependence effects on the population at each timestep
#' @param dynamics_order the order in which the population dynamics should be executed
#'  on the landscape object - default is "change" -> "dispersal" -> "modification" -> "density_dependence". 
#'  Note, if population dynamics are reordered, all dynamics must be listed in \code{dynamics_order}.
#'
#'
#' @return An object of class \code{population_dynamics}
#' 
#' @export
#'
#' @examples
#' 
#' # Example of setting up population dynamics to only use a population change function.
#' 
#' \dontrun{
#' ls <- landscape(population = egk_pop, suitability = NULL, carrying_capacity = NULL)
#' 
#' pd <- population_dynamics(change = growth(egk_mat))
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

population_dynamics <- function (change = NULL,
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL,
                                 dynamics_order = c("change", "dispersal", "modification", "density_dependence")) {
  
  if (!is.null(density_dependence)) {
    
    # transfer information about the density dependence- contributing stages to the functions that need it
    dd_stages <- attr(density_dependence, "density_dependence_stages")
    
    if (!is.null(dispersal)) {
      # create an object called "density_dependence_stages" in the environment
      # in which dispersal runs
      environment(dispersal)$density_dependence_stages <- dd_stages
    }
    
  }
  
  order_sorted <- sort(dynamics_order)
  order_sorted_target <- c("change", "density_dependence", "dispersal", "modification")
  
  if (!identical(order_sorted, order_sorted_target)) {
    stop ("You must specify all population dynamics regardless of whether they are\n",
          "used or not (set to NULL). Please correct and re-run the simulation.")
  }
  
  pop_dynamics <- function (landscape, timestep) {

    for (i in dynamics_order) {
      if (is.null(get(i))) next
      landscape <- do.call(i, list(landscape, timestep))
    }
      
    landscape
  }
  
  as.population_dynamics(pop_dynamics)
  
}

is.population_dynamics <- function (x) {
  inherits(x, 'population_dynamics')
}

# #' @rdname population_dynamics
# #'
# #' @export
# #'
# #' @examples
# #'
# #' # Print details about the 'population_dynamics' object 
# #' print(pd)
# 
# print.population_dynamics <- function (x, ...) {
#   cat("This is a population_dynamics object.")
# }

##########################
### internal functions ###
##########################

as.population_dynamics <- function (population_dynamics_function) {
  as_class(population_dynamics_function, "population_dynamics", "function")
}
