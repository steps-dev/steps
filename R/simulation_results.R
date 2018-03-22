#' Run an simulation experiment to make many spatially-explicit population projections and look at the variation.
#'
#' @description A habitat object is used to store spatially-explicit information on habitat suitability and the carrying_capacity of a landscape.
#' It is a sub-component of a \link[steps]{state} object and is modified in each timestep of an experiment.
#' 
#' @rdname simulation_results
#'
#' @param state a state object - static habitat, population, and demography in a timestep
#' @param dynamics a dynamics object - modules that change habitat, population, and demography during and experiment
#' @param timesteps number of timesteps used in the experiment
#' @param simulations an experiment_reults object
#' 
#' @return An object of class \code{simulation_results}
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
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"), carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(transition_matrix = steps:::fake_transition_matrix(4), dispersal_parameters = rlnorm(1))
#' test_population <- build_population(stack(replicate(4, test_habitat$carrying_capacity * 0.2)))
#' test_state <- build_state(test_habitat, test_demography, test_population)
#' fast_approximation <- build_dynamics(steps:::no_habitat_dynamics, steps:::no_demographic_dynamics, steps:::fast_population_dynamics)
#' simulation_results <- simulation(test_state,fast_approximation,timesteps=10,simulations=10)
 
simulation <- function(state, dynamics, timesteps, simulations, check=TRUE, ...){
                         
                          if(check) {
                            message('running a single test to see if the experiment works, simulation will stop if this fails\n')
                            stopifnot(is.experiment_results(experiment(state,dynamics,timesteps=1)))
                          }
                          
                          future::plan(multiprocess)
                          simulation_results <- list()
                          for(sim in seq_len(simulations)){
                            simulation_results[[sim]] <- future(
                              experiment(state,dynamics,timesteps)
                            )
                          }
                          
                          set_class(simulation_results, "simulation_results")
}

#' @rdname simulation_results
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'simulation results'
#'   
#' is.simulation_results(results)

is.simulation_results <- function (results) {
  inherits(results, 'simulation_results')
}

