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
#' library(steps)
#' library(raster)
#' 
#' # Import a raster layer for habitat
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' # Create a life-stage matrix
#' mat <- matrix(c(0.000,0.000,0.302,0.302,
#'                 0.940,0.000,0.000,0.000,
#'                 0.000,0.884,0.000,0.000,
#'                 0.000,0.000,0.793,0.793),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#' 
#' # Create a matrix with standard deviations for environmental stochasticity
#' mat_sd <- matrix(c(0.000,0.00,0.010,0.010,
#'                 0.010,0.000,0.000,0.000,
#'                 0.000,0.010,0.000,0.000,
#'                 0.000,0.000,0.010,0.010),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat_sd) <- rownames(mat_sd) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#' 
#' # Create a stack of raster layers to represent each
#' # life-stage of a population structure (four in this case)
#' pop <- stack(replicate(4, ceiling(r * 0.2)))
#'
#' # Create raster and shuffle values (omit NAs)
#' r2 <- r
#' r2[na.omit(r2)] <- sample(r[na.omit(r)])
#' 
#' # Create raster and shuffle values (omit NAs)
#' r3 <- r
#' r3[na.omit(r3)] <- sample(r[na.omit(r)])
#' 
#' # Create a list of rasters stacks for all life stages
#' surv <- list(stack(r2, r2, r2),
#'              stack(r2, r2, r2),
#'              stack(r2, r2, r2),
#'              stack(r2, r2, r2))
#'
#' # Create a list of raster stacks when the first two stages are NULL            
#' fec <- list(NULL,
#'             NULL,
#'             stack(r3, r3, r3),
#'             stack(r3, r3, r3))
#' 
#' # Construct habitat, demography, and population objects.
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"),
#'                               carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(transition_matrix = mat)
#' test_population <- build_population(pop)
#' 
#' # Construct a state object
#' test_state <- build_state(test_habitat, test_demography, test_population)
#'
#' # Select existing habitat functions to be run on the demography
#' # in a simulation and specify input parameters: 
#' fires <- disturbance_fires(habitat_suitability = r / cellStats(r, "max"),
#'                            disturbance_layers = dist,
#'                            effect_time = 1)
#' 
#' # Construct a habitat dynamics object
#' test_hab_dynamics <- build_habitat_dynamics(fires)

build_habitat_dynamics <- function (...) {
  
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