#' How the population changes in a landscape.
#'
#' Pre-defined or custom functions to define population change during a simulation.
#' Please see the tutorial vignette titled "Creating custom *steps* functions"
#' for information on how to write custom functions for use in simulations.
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
#' throughout a simulation. Please see the tutorial vignette titled "Creating custom
#' *steps* functions" for information on how to write custom functions for use in simulations.
#'
#' @param transition_matrix A symmetrical age-based (Leslie) or stage-based (Lefkovitch)
#'   population structure matrix.
#' @param global_stochasticity,local_stochasticity Either scalar values or
#'   matrices (with the same dimension as \code{transition_matrix}) specifying
#'   the variability (in standard deviations) in the transition matrix either for
#'   populations in all grid cells (\code{global_stochasticity}) or for each
#'   grid cell population separately (\code{local_stochasticity})
#' @param transition_function A function to specify or modify life-stage transitions
#'   at each timestep. See \link[steps]{transition_function}.
#' @param transition_order Order of transitions performed in growth function. This behaviour
#'   is only applied when demographic stochasticity is set to "full" (default) and transitions
#'   are applied sequentially. By default "fecundity" is performed first (calculating the
#'   number of new individuals to be added to the populations), then "survival" is applied.
#'   The final population is the sum of these. Users should be cautious of specifying
#'   "survival" to be performed first as typically survival of reproductive stages will already
#'   be accounted for in the fecundity values of the transition matrix.
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
                    transition_function = NULL,
                    transition_order = c("fecundity", "survival")) {
  
  transition_order <- match.arg(transition_order)
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
    
    # calculate global and local noise
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
    
    values <- pmax_zero(values)
    values <- pmin(values, rep(upper, n_cells))
    transition_array[idx_full] <- values
    
    if (steps_stash$demo_stochasticity == "full") {
      
      total_pop <- rowSums(population)
      has_pop <- total_pop > 0

      # browser()
      
      if (transition_order == "fecundity" && sum(has_pop) >= 1) {
        # first step - perform fecundity to add individuals to the populations
        new_population <- add_offspring(population[has_pop, ], transition_array[ , , has_pop])
        # second step - perform survival on new population
        surv_population <- surviving_population(population[has_pop, ], transition_array[ , , has_pop])
        # add new and surviving populations
        surv_population[ , 1] <- surv_population[ , 1] + new_population
        population[has_pop, ] <- surv_population
      }
      
      if (transition_order == "survival" && sum(has_pop) >= 1) {
        # first step - perform survival on population
        surv_population <- surviving_population(population[has_pop, ], transition_array[ , , has_pop])
        # second step - perform fecundity to add individuals to the populations
        new_population <- add_offspring(surv_population, transition_array[ , , has_pop])
        # add new and surviving populations
        surv_population[ , 1] <- surv_population[ , 1] + new_population
        population[has_pop, ] <- surv_population
      }
      
      
    } else {
      
      population <- sapply(seq_len(n_cells),
                           function(x) transition_array[ , , x] %*% matrix(population[x, ]))
      
      population <- t(population)

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

add_offspring <- function (population, transition_array) {

  pops <- population
  
  # get fecundities for all eligible stages
  # N.B. assumes we can only recruit into the first stage!
  if (class(transition_array) == "matrix") {
    fecundities <- t(transition_array[1, ])
  } else {
    fecundities <- t(transition_array[1, , ])
  }

  
  # get expected number, then do a poisson draw about this
  expected_offspring <- fecundities * pops
  new_offspring_stochastic <- expected_offspring
  new_offspring_stochastic[] <- stats::rpois(length(expected_offspring), expected_offspring[])
  
  # sum stage 1s created by all other stages
  new_offspring <- rowSums(new_offspring_stochastic)
  return(new_offspring)
}

surviving_population <- function (population, transition_array) {
  survival_array <- transition_array
  
  if (class(transition_array) == "matrix") {
    survival_array[1, ] <- 0
  } else {
    survival_array[1, , ] <- 0
  }

  if(inherits(population, "numeric")) {
    n_stage <- length(population)
  } else {
    n_stage <- ncol(population)
  }
   
  
  # loop through stages, getting the stages to which they move (if they survive)
  if(inherits(population, "numeric")) {
    survival_stochastic <- matrix(0, 1, n_stage)
  } else {
    survival_stochastic <- matrix(0, nrow(population), n_stage)
  }

  for (stage in seq_len(n_stage)) {
    
    # get the populations that have any individuals of this stage
    if(inherits(population, "numeric")) {
      pops <- population[stage]
    } else {
      pops <- population[, stage]
    }
    
    # probability of transitioning to each other stage
    if (class(survival_array) == "matrix") {
      probs <- t(survival_array[ , stage])
    } else {
      probs <- t(survival_array[ , stage, ])
    }
    
    # add on probability of dying
    
    surv_prob <- rowSums(probs)
    
    # check for sensible values
    if(any(surv_prob > 1)) {
      stop("Survival values greater than one have been detected. Please check to ensure\n",
           "that your transition matrix values or stochasticity (standard deviation) values\n",
           "are sensible and re-run the simulation.")
    }
    
    probs <- cbind(probs, 1 - surv_prob)
    
    # loop through cells with population (rmultinom is not vectorised on probabilities)
    new_stages <- matrix(NA, length(pops), n_stage)
    idx <- seq_len(n_stage)
    for (i in seq_len(length(pops))) {
      new_stages[i, ] <- stats::rmultinom(1, pops[i], probs[i, ])[idx, ]
    }
    
    # update the population
    survival_stochastic <- survival_stochastic + new_stages
    
  }
  
  return(survival_stochastic)
}
