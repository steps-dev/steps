#' How the population changes in a landscape.
#'
#' Pre-defined or custom functions to define population change during a simulation.
#'
#' @name population_change_functions
#'
#' @seealso
#' \itemize{
#'   \item{\link[steps]{growth} is a default function for changing populations based on
#'   transition matrices and functions}
#'   }
NULL

#' Population growth
#' 
#' This function applies negative or positive growth to the population using matrix
#' multiplication. Stochasticty can be added to cell-based transition matrices or globally.
#' Users can also specify a built-in or custom function to modify the transition matrices
#' throughout a simulation.
#'
#' @param transition_matrix A symmetrical age-based (Leslie) or stage-based (Lefkovitch)
#'   population structure matrix.
#' @param global_stochasticity,local_stochasticity either scalar values or
#'   matrices (with the same dimension as \code{transition_matrix}) specifying
#'   the variability (in standard deviations) in the transition matrix either for
#'   populations in all grid cells (\code{global_stochasticity}) or for each
#'   grid cell population separately (\code{local_stochasticity})
#' @param transition_function A custom function defined by the user specifying
#'   modifications to life-stage transitions at each timestep. See \link[steps]{transition_function}.
#' 
#' @export
#' 
#' @examples
#' 
#' # Example of a growth function that changes the populations based on a transition matrix that
#' # is subject to global stochasticity. 
#' 
#' \dontrun{
#' stoch_growth <- growth(transition_matrix = egk_mat, global_stochasticity = egk_mat_stoch)
#' 
#' ls <- landscape(population = egk_pop, suitability = NULL, carrying_capacity = NULL)
#' 
#' pd <- population_dynamics(change = stoch_growth)
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

growth <- function (transition_matrix,
                    global_stochasticity = 0,
                    local_stochasticity = 0,
                    transition_function = NULL) {
  
  is_function <- inherits(transition_function, "function")
  
  # if it's a list of functions, chain them together
  if (!is_function && inherits(transition_function, "list")) {
    
    # check the elements of this list are all functions, with the right arguments
    all_funs <- all(unlist(lapply(transition_function,
                                 function(x) inherits(x, "function"))))
    expected_params <- all(unlist(lapply(transition_function,
                                         function(x) names(formals(x)) %in% c("transition_array",
                                                                              "landscape",
                                                                              "timestep"))))
    
    if (!all_funs || !expected_params) {
      stop("A transition function list must contain only function objects that\n",
           "each include the arguments: transition_array, landscape, and timestep.\n",
           "Please check your inputs and re-run the simulation.")
    }
    
    transition_function_list <- transition_function
    
    
    transition_function <- function(transition_array, landscape, timestep) {

      for (fun in transition_function_list) {
        transition_array <- fun(transition_array, landscape, timestep)
      }
      
      transition_array
      
    }
    
    is_function <- TRUE
    
  }
  
  idx <- which(transition_matrix != 0)
  is_recruitment <- upper.tri(transition_matrix)[idx]
  upper <- ifelse(is_recruitment, Inf, 1)
  vals <- transition_matrix[idx]
  dim <- nrow(transition_matrix)
  
  if (is.matrix(global_stochasticity)) {
    stopifnot(identical(c(dim, dim), dim(global_stochasticity)))
    stopifnot(identical(which(global_stochasticity != 0), idx))
    global_stochasticity <- global_stochasticity[idx]
  }

  if (is.matrix(local_stochasticity)) {
    stopifnot(identical(c(dim, dim), dim(local_stochasticity)))
    stopifnot(identical(which(local_stochasticity != 0), idx))
    local_stochasticity <- local_stochasticity[idx]
  }
  
  pop_dynamics <- function (landscape, timestep) {
    
    # import components from landscape object
    population_raster <- landscape$population
    
    # get population as a matrix
    cell_idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population <- raster::extract(population_raster, cell_idx)
    
    n_cells <- length(cell_idx)
    
    # add global noise to the transition matrices and truncate
    global_noise <- stats::rnorm(length(idx), 0, global_stochasticity)
    local_noise <- stats::rnorm(length(idx) * n_cells, 0, local_stochasticity)
    total_noise <- global_noise + local_noise

    # pad the index to get corresponding elements in each slice  
    addition <- length(transition_matrix) * (seq_len(n_cells) - 1)
    idx_full <- as.numeric(outer(idx, addition, FUN = "+"))
    
    if (is_function) {
      
      # create transition array and fill with initial matrix values
      transition_array <- array(0, c(dim, dim, n_cells))
      transition_array[] <- transition_matrix[]
      
      # update the transition array
      transition_array <- transition_function(transition_array, landscape, timestep)
      
      values <- transition_array[idx_full] + total_noise
    } else {
      transition_array <- array(0, c(dim, dim, n_cells))
      values <- vals + total_noise
    }
    
    values <- pmax(values, 0)
    values <- pmin(values, rep(upper, n_cells))
    transition_array[idx_full] <- values

    if (steps_stash$demo_stochasticity == "full") {

      n_cell <- nrow(population)
      n_stage <- ncol(population)
      
      survival_array <- transition_array
      survival_array[1, , ] <- 0

      # loop through stages, getting the stages to which they move (if they survive)
      survival_stochastic <- matrix(0, n_cell, n_stage)
      for (stage in seq_len(n_stage)) {
        
        # get the populations that have any of this stage
        pop <- population[, stage]
        any <- pop > 0
        n_any <- sum(any)
        
        # how many to transition
        sizes <- pop[any]

        # probability of transitioning to each other stage
        probs <- t(survival_array[, stage, any])

        # add on probability of dying
        surv_prob <- rowSums(probs)
        probs <- cbind(probs, 1 - surv_prob)

        # loop through cells with population (rmultinom is not vectorised on probabilities)
        new_stages <- matrix(NA, n_any, n_stage)
        idx <- seq_len(n_stage)
        for (i in seq_len(n_any)) {
          new_stages[i, ] <- stats::rmultinom(1, sizes[i], probs[i, ])[idx, ]
        }
        
        # update the population
        survival_stochastic[any, ] <- survival_stochastic[any, ] + new_stages

      }
      
      # do reproduction all in one go
      
      # fecundity
      
      # only look at sites with any individuals
      total_pop <- rowSums(population)
      any <- total_pop > 0
      n_any <- sum(any)

      pop_any <- population[any, ]
      
      # get fecundities
      # N.B. assumes we can only recruit into the first stage!
      fecundities <- t(transition_array[1, , any])

      # get expected number, then do a poisson draw about this
      expected_offspring <- fecundities * pop_any
      new_offspring_stochastic <- expected_offspring
      new_offspring_stochastic[] <- stats::rpois(length(expected_offspring), expected_offspring[])

      # sum stage 1s created by all other stages
      new_offspring <- rowSums(new_offspring_stochastic)

      # combinethe survivals and fecundities
      population <- survival_stochastic
      population[any, 1] <- population[any, 1] + new_offspring

    } else {
      
      population <- sapply(seq_len(n_cells),
                           function(x) transition_array[ , , x] %*% matrix(population[x, ]))
      
      population <- t(population)
      
      # get whole integers
      population <- round_pop(population)

    }
    
    # put back in the raster
    population_raster[cell_idx] <- population
    
    landscape$population <- population_raster
    
    landscape
  }
  
  result <- as.population_growth(pop_dynamics)

  result
}


##########################
### internal functions ###
##########################

as.population_growth <- function (simple_growth) {
  as_class(simple_growth, "population_growth", "function")
}