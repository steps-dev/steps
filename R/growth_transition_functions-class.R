#' Create a growth transition function
#'
#' A growth transition function defines how spatial objects or custom functions influence
#' survival and fecundity. Two built-in functions are provided for the user to select, however,
#' a user may also provide custom written functions to modify survival and fecundity throughout
#' a simulation. Please see the tutorial vignette titled "Creating custom *steps* functions"
#' for information on how to write custom functions for use in simulations.
#' 
#' 
#' @name transition_function
#' @seealso
#' \itemize{
#'   \item{\link[steps]{modified_transition} to use rasters to modify survival and fecundity}
#'   \item{\link[steps]{competition_density} to use relationship to carrying capacity to modify
#'   survival and fecundity}
#'   }
NULL

#' Spatially-explicit transition function
#' 
#' In the built-in \code{modified_transition function}, the values of fecundity and survival
#' in local cell-based transition matrices are multiplied by values in the named spatial objects
#' for each cell. The spatial objects can be rasters that are stored in the landscape object.
#' 
#' The behaviour of the function is to modify any non-zero values in the first row by
#' the "fecundity_layer" and non-zero values in rows other than the first by the "survival_layer".
#' This is irrespective of the type of matrix or any assumptions made by the user in creating
#' the transition matrix. For example, if the transition matrix values include both the
#' probabilities of surviving AND growing into the next stage, these can NOT be modified
#' individually. This operation would require the use of a custom function - see the "Creating
#' custom *steps* functions" vignette for more information.
#'
#' @param survival_layer the name of a spatial layer in the landscape object used to modify survival values (i.e. non-zero values in the first row).
#' @param fecundity_layer the name of a spatial layer in the landscape object used to modify fecundity values (i.e. non-zero values in rows other than the first).
#' 
#' @return An object of class \code{transition_function}
#' 
#' @export
#'
#' @examples
#' 
#' # Vital rates (survival and fecundity) modified based on habitat suitability.
#' 
#' \dontrun{
#' mod_fun <- modified_transition(survival_layer = "suitability", fecundity_layer = "suitability")
#' 
#' ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = NULL)
#' 
#' pd <- population_dynamics(change = growth(egk_mat, transition_function = mod_fun))
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

modified_transition <- function(survival_layer = NULL,
                                fecundity_layer = NULL) {
  
  fun <- function (transition_array, landscape, timestep) {
    
    transition_matrix <- transition_array[, , 1]
    idx <- which(transition_matrix != 0)
    is_recruitment <- upper.tri(transition_matrix)[idx]
    
    array_length <- dim(transition_array)[3]

    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
    
    if (is.null(survival_layer)) {
      surv_mult <- rep(1, length(cell_idx))
    } else {
      if (raster::nlayers(landscape$suitability) > 1) {
        surv_mult <- landscape[[survival_layer]][[timestep]][cell_idx]
      } else {
        surv_mult <- landscape[[survival_layer]][cell_idx]
      }
    }
    
    if (is.null(fecundity_layer)) {
      fec_mult <- rep(1, length(cell_idx))
    } else {
      if (raster::nlayers(landscape$suitability) > 1) {
        fec_mult <- landscape[[fecundity_layer]][[timestep]][cell_idx]
      } else {
        fec_mult <- landscape[[fecundity_layer]][cell_idx]
      }
    }
    
    for (i in seq_len(array_length)) {
      transition_array[, , i][idx[!is_recruitment]] <- transition_array[, , i][idx[!is_recruitment]] * surv_mult[i]
      transition_array[, , i][idx[is_recruitment]] <- transition_array[, , i][idx[is_recruitment]] * fec_mult[i]
    }

    transition_array
    
  }
  
  as.transition_function(fun)
  
}

#' Competition density function
#'
#' Adjusts the life-stage transition matrix in each cell based on the carrying capacity in the cell and
#' a density dependence function - default is Beverton-Holt. The user may specify which life-stages are 
#' affected by density dependence. If \code{R_max} is not provided this is calculated from the local cell-based
#' transition matrices internally. By providing initial stable age distribution values, performance can be
#' increased as the function internally calculates these values through optimisation.
#' 
#' @param stages which life-stages contribute to density dependence - default is all
#' @param mask a matrix of boolean values (TRUE/FALSE), equal in dimensions to the life-stage transition matrix
#' and specifying which vital rates (i.e. survival and fecundity) are to be modified by the function
#' @param R_max optional value of maximum growth rate (lambda) if known
#' @param stable_age optional vector of stable age distributions if known
#' 
#' @export
#'
#' @examples
#' 
#' # Vital rates (survival and fecundity) modified based on approach to carrying capacity
#' # by the 2nd and 3rd life stages.
#' 
#' \dontrun{
#' mod_fun <- competition_density(stages = c(2, 3))
#' 
#' ls <- landscape(population = egk_pop, suitability = NULL, carrying_capacity = egk_k)
#' 
#' pd <- population_dynamics(change = growth(egk_mat, transition_function = mod_fun))
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

competition_density <- function(stages = NULL,
                                mask = NULL,
                                R_max = NULL,
                                stable_age = NULL) {

  fun <- function (transition_array, landscape, timestep) {

    # get metrics and constructor info
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
    n_cells <- length(cell_idx)
    
    # get population matrix
    pop_raster <- landscape$population
    population <- raster::extract(pop_raster, cell_idx)
    
    # get carrying capacity (internal function to STEPS)
    # 22.01.20 - # cc <- get_carrying_capacity(landscape, timestep)
    # 22.01.20 - # K <- raster::extract(cc, cell_idx)
    K <- raster::extract(landscape$carrying_capacity, cell_idx) # 22.01.20
    
    if (!is.null(stages)) {
      if (length(stages) == 1) {
        N <- population[, stages]
      } else {
        N <- rowSums(population[, stages])
      }
    } else {
      N <- rowSums(population)
    }

    target_cells <- which(N - K != 0 & N != 0)

    # modify life-stage transition array
    for (i in target_cells) {
      transition_array[, , i] <- density_modified_transition(N = N[i],
                                                             K = K[i],
                                                             transition_matrix = transition_array[, , i],
                                                             mask = mask,
                                                             R_max = R_max,
                                                             stable_age = stable_age)
    }
    
    # return array with required dimensions
    transition_array
    
  }
  
  as.transition_function(fun)
  
}


# #' @rdname transition_function
# #'
# #' @param x an object to print or test as a transition_function object
# #' @param ... further arguments passed to or from other methods
# #'
# #' @export
# #'
# #' @examples
# #'
# #' print(test_transition_function)
# 
# print.transition_function <- function (x, ...) {
#   cat("This is a transition_function object")
# }

##########################
### internal functions ###
##########################

as.transition_function <- function (transition_function) {
  as_class(transition_function, "transition_function", "function")
}

get_R <- function (transition_matrix, stable_age = NULL, tolerance = 0.001, max_iter = 100) {
  
  if (is.null(stable_age)) {
    stable_age <- rep(1, ncol(transition_matrix))
  }
  
  old_stages <- stable_age
  converged <- FALSE
  iter <- 0
  old_Rs <- rep(.Machine$double.eps, ncol(transition_matrix))

  while (!converged & iter < max_iter) {
    new_stages <- transition_matrix %*% old_stages
    Rs <- new_stages / old_stages
    errors <- abs(1 - (Rs / old_Rs))
    converged <- all(errors < tolerance)
    old_Rs <- Rs
    old_stages <- new_stages
    iter <- iter + 1
  }
  
  warn_once(!converged,
            paste("estimation of growth rate did not converge in",
                  max_iter,
                  "iterations"),
            warning_name = "growth_estimation")

  # return the intrinsic growth rate
  Rs[1]
  
}

ideal_R <- function (K, N, R_max) {
  if(R_max > 1) {
    num <- R_max * K
    denom <- R_max * N - N + K
    R <- num / denom
  } else {
    R <- R_max
  }
  R
}

# multiply m by the relevant elements of transition_matrix (specified by mask) and return the growth rate R
apply_m <- function (m, transition_matrix, mask = NULL) {
  if (is.null(mask)) {
    transition_matrix <- transition_matrix * m
  } else {
    mask <- as.logical(mask)
    transition_matrix[mask] <- transition_matrix[mask] * m
  }
  transition_matrix
}

# find a value of m with which to modify transition_matrix, to get to this target value of R
find_m <- function(R_target, transition_matrix, mask = NULL, stable_age = NULL) {
  
  obj <- function (m, R_target, transition_matrix, mask = NULL, stable_age = NULL) {
    new_transition_matrix <- apply_m(m, transition_matrix, mask)
    R_current <- get_R(new_transition_matrix, stable_age = stable_age)
    (R_current - R_target) ^ 2
  } 
  
  out <- stats::optimise(f = obj,
                         interval = c(0, 5),
                         R_target,
                         transition_matrix,
                         mask,
                         stable_age)
  out$minimum
  
}

density_modified_transition <- function (N,
                                         K,
                                         transition_matrix,
                                         R_max = NULL,
                                         stable_age = NULL,
                                         mask = NULL) {
  
  # if the optimal R isn't provided, recalculate it (ideally pre-calculate it to
  # save computation)
  init_Rmax_null <- is.null(R_max)
  
  if (init_Rmax_null) {
    R_max <- get_R(transition_matrix, stable_age = stable_age)
  }
  
  # get the target value of R, for this degree of over/under-population
  if (R_max > 1) {
    R_target <- ideal_R(K, N, R_max)
    #print(R_target)
    
    # find a value of m with which to modify transition_matrix, to get to this target value of R
    m <- find_m(R_target, transition_matrix, mask, stable_age = stable_age)
    #print(m)
    
    # multiply m by the relevant bits of transition_matrix
    transition_matrix <- apply_m(m, transition_matrix, mask)
    
  }
  
  transition_matrix
  
}
