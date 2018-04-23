#' Change the demography in a state object
#' 
#' @description A 'demography dynamics' object is used to modify life-stage transition matrices - adding stochasticity for example.
#' It is a sub-component of a \link[steps]{dynamics} object and is executed in each timestep of an experiment.
#'
#' @rdname demography_dynamics
#'
#' @param demography_dynamics_function A function that operates on a state object to change demography at specified timesteps. User may enter a custom function or select a pre-defined module - see examples. 
#' @param x an object to print or test as an demography_dynamic object
#' @param ... further arguments passed to or from other methods
#' @param env_stoch a function for adding environmental stochasticity to the global transition matrix at each timestep
#' @param demo_dens_dep a function for modifying the transition matrix at each timestep when carrying capacity is reached
#' @param global_transition_matrix a life-stage transition matrix
#' @param stochasticity a matrix with standard deviations (consistent or varying) around the transition means with matching dimensions as the life-stage transition matrix or a number representing a consitent standard deviation to apply to all transitions (default is 0)
#' @param fecundity_fraction a multiplier value between 0 and 1 for fecundity values in the transitiopn matrix
#' @param survival_fraction a multiplier value between 0 and 1 for fecundity values in the transitiopn matrix 
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
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' mat <- matrix(c(0.000,0.000,0.302,0.302,
#'                 0.940,0.000,0.000,0.000,
#'                 0.000,0.884,0.000,0.000,
#'                 0.000,0.000,0.793,0.793),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#' 
#' mat_sd <- matrix(c(0.000,0.000,1,1,
#'                 1,0.000,0.000,0.000,
#'                 0.000,1,0.000,0.000,
#'                 0.000,0.000,1,1),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat_sd) <- rownames(mat_sd) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#' 
#' pop <- stack(replicate(4, ceiling(r * 0.2)))
#' 
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"),
#'                               carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(transition_matrix = mat,
#'                                     dispersal_parameters = rlnorm(1))
#' test_population <- build_population(pop)
#' 
#' test_state <- build_state(test_habitat, test_demography, test_population)
#' 
#' # Create a generic function that simply returns an unmodified state
#' # object at each timestep
#'
#' example_function <- function (state, timestep) {
#'   state
#' }
#' 
#' # Define the function as a demography_dynamics object
#' 
#' example_function <- as.demography_dynamics(example_function)

as.demography_dynamics <- function (demography_dynamics_function) {
  as_class(demography_dynamics_function, "demography_dynamics", "function")
}

#' @rdname demography_dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'demography dynamics'
#'   
#' is.demography_dynamics(example_function)

is.demography_dynamics <- function (x) {
  inherits(x, 'demography_dynamics')
}

#' @rdname demography_dynamics
#'
#' @export
#'
#' @examples
#' 
#' # Print details about the demography_dynamics object
#' 
#' print(example_function)

print.demography_dynamics <- function (x, ...) {
  cat("This is a demography_dynamics object")
}


#' @rdname demography_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the demography_dynamics function to modify the transition
#' # matrix:
#' 
#' env_stoch <- demo_environmental_stochasticity(global_transition_matrix = mat,
#'                                              stochasticity = mat_sd)
#' 
#' example_function <- demography_dynamics(env_stoch,
#'                                         demo_dens_dep =  demo_density_dependence())

demography_dynamics <- function (env_stoch = NULL, demo_dens_dep = NULL) {
  
  demographic_dynamics <- function (state, timestep) {
    
    if (!is.null(env_stoch))
      state <- env_stoch(state, timestep)
    
    if (!is.null(demo_dens_dep))
      state <- demo_dens_dep(state, timestep)

    state
    
  }
  
  as.demography_dynamics(demographic_dynamics)
  
}

##########################
### internal functions ###
##########################

as.demography_environmental_stochasticity <- function (demography_environmental_stochasticity) {
  as_class(demography_environmental_stochasticity, "demography_dynamics", "function")
}

as.demography_density_dependence <- function (demography_density_dependence) {
  as_class(demography_density_dependence, "demography_dynamics", "function")
}


####################################
### pre-defined module functions ###
####################################

#' @rdname demography_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the demo_environmental_stochasticity function to modify the transition
#' # matrix with specified environmental stochasticity:
#' 
#' test_demo_es <- demo_environmental_stochasticity(global_transition_matrix = mat,
#'                                     stochasticity = mat_sd)

demo_environmental_stochasticity <- function (global_transition_matrix,
                                              stochasticity=0) {
  
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
  
  env_stoch_fun <- function (state, timestep) {
    
    transition_matrix <- global_transition_matrix
    
    transition_matrix[idx] <- extraDistr::rtnorm(length(idx),
                                                 vals,
                                                 stochasticity,
                                                 a = lower,
                                                 b = upper)
    
    state$demography$transition_matrix <- transition_matrix
    
    state
    
  }
  
  as.demography_environmental_stochasticity(env_stoch_fun)
  
}


#' @rdname demography_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the demo_density_dependence function to modify the transition
#' # matrix once carrying capacity is reached:
#' 
#' test_demo_dd <- demo_density_dependence(fecundity_fraction = 1,
#'                                                 survival_fraction = 0.5)

demo_density_dependence <- function (fecundity_fraction = 1,
                               survival_fraction = 1) {
  
  dens_dep_fun <- function (state, timestep) {
    
    idx <- which(!is.na(raster::getValues(state$population$population_raster[[1]])))
    population <- raster::extract(state$population$population_raster, idx)
    
    transition_matrix <- state$demography$transition_matrix
    
    if (any(rowSums(population) > as.vector(state$habitat$carrying_capacity))) {

      ids <- which(transition_matrix != 0 & row(transition_matrix) == 1)
      transition_matrix[ids] <- transition_matrix[ids] * fecundity_fraction
      
      
      ids <- which(transition_matrix != 0 & row(transition_matrix) != 1)
      transition_matrix[ids] <- transition_matrix[ids] * survival_fraction

      state$demography$transition_matrix <- transition_matrix
      
    }
    
    state
    
  }
  
  as.demography_density_dependence(dens_dep_fun)
  
}


