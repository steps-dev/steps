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
#' @param state a state object to apply the demographic function to
#' @param timestep the timestep in the experiment to apply the demographic function to the state object
#' @param global_transition_matrix a life-stage transition matrix
#' @param stochasticity a matrix with standard deviations (consistent or varying) around the transition means with matching dimensions as the life-stage transition matrix or a number representing a consitent standard deviation to apply to all transitions (default is 0)
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
#' state
#' }
#' 
#' # Define the function as a demography_dynamics object
#' 
#' no_demography_dynamics <- as.demography_dynamics(example_function)
#' 
#' # Alternatively, embedded in the function:
#' 
#' example_function <- function(...) {
#'   int.func <- function (state, timestep) {
#'     state
#'    }
#' as.demography_dynamics(int.func)
#' }

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
#' # Print details about the demography_dynamics object
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

#' @rdname demography_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the no_demography_dynamics object as a placeholder as it 
#' # does not modify the demography object:
#' 
#' test_state2 <- no_demography_dynamics(test_state, 1)
#' 
#' identical(test_state, test_state2)
 
no_demography_dynamics <- function () {

    demographic_dynamics <- function (state, timestep) {
    state
    }

  as.demography_dynamics(demographic_dynamics)

}


#' @rdname demography_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the envstoch_demography_dynamics function to modify the transition
#' # matrix with specified environmental stochasticity:
#' 
#' test_state2 <- envstoch_demography_dynamics(global_transition_matrix = mat,
#'                                     stochasticity = mat_sd)

envstoch_demography_dynamics <- function (global_transition_matrix, stochasticity=0) {
  
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