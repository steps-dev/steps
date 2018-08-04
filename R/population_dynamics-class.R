#' Modify the population in a state object.
#' 
#' A \code{population_dynamics} object is used to modify species populations
#' in space and time.
#' 
#' A \code{population_dynamics} object is a sub-component of a \link[steps]{dynamics}
#' object and is executed in each timestep of a simulation.  Note, some dynamics
#' functions can be executed at non-regular intervals (i.e. only timesteps
#' explicitly defined by the user). The \code{build_population_dynamics} function is
#' used to construct a population dynamics object consisting of several population
#' dynamics functions and their associated parameters. These functions specify how
#' the population in the state object will be modified throughout a simulation. 
#'
#' @rdname population_dynamics
#'
#' @param pop_change a function to define how population growth occurs at each timestep
#' @param pop_disp a function to define how the population disperses at each timestep (default is exponential kernel)
#' @param pop_mod a function to define any deterministic changes to the population - such as translocation - at each timestep
#' @param pop_dens_dep a function to control density dependence effects on the population at each timestep
#' @param object a population_dynamic object
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

#' # Use the population_dynamics object to modify the population with
#' # a default population change function:
#' test_pop_dynamics <- build_population_dynamics()
#' test_state2 <- test_pop_dynamics(test_state, 1)
#' 
#' par(mfrow=c(1,2))
#' plot(test_state$population$population_raster[[2]])
#' plot(test_state2$population$population_raster[[2]])

build_population_dynamics <- function (pop_change = simple_growth(),
                                 pop_disp = NULL,
                                 pop_mod = NULL,
                                 pop_dens_dep = NULL) {
  
  pop_dynamics <- function (state, timestep) {
    
    state <- pop_change(state, timestep)
    
    if (!is.null(pop_disp))
      state <- pop_disp(state, timestep)
    
    if (!is.null(pop_mod))
      state <- pop_mod(state, timestep)
    
    if (!is.null(pop_dens_dep))
      state <- pop_dens_dep(state, timestep)
    
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

is.population_dynamics <- function (object) {
  inherits(object, 'population_dynamics')
}

#' @rdname population_dynamics
#'
#' @export
#'
#' @examples
#'
#' # Print details about the 'population_dynamics' object 
#' print(test_pop_dynamics)

print.population_dynamics <- function (object, ...) {
  cat("This is a population_dynamics object")
}