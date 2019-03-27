#' Create a growth transition function
#'
#' @description A growth transition function defines how spatial objects or custom functions influence
#' survival and fecundity. A user may select from built-in functions or provide a custom written function
#' to modify survival and fecundity throughout a simulation.
#' 
#' In the built-in \code{modified_transition function}, the values of fecundity and survival
#' in local cell-based transition matrices are multiplied by values in the named spatial objects
#' for each cell. The spatial objects can be rasters that are stored in the landscape object.
#' 
#' A commonly used \code{competition_density} dependence function is also provided in the software for
#' the user to select, however, a user may also provide other custom written density dependence functions.
#' 
#' @rdname transition_function
#'
#' @param transition_matrix a symmetrical age-based (Leslie) or stage-based population
#'   structure matrix.
#' @param survival_layer the name of a spatial layer in the landscape object used to modify survival values.
#' @param fecundity_layer the name of a spatial layer in the landscape object used to modify fecundity values.
#' 
#' @return An object of class \code{transition_function}
#' 
#' @export
#'
#' @examples
#' 
#' test_mod_transition <- modified_transition(egk_mat)

modified_transition <- function(transition_matrix,
                                survival_layer = NULL,
                                fecundity_layer = NULL) {
  
  idx <- which(transition_matrix != 0)
  is_recruitment <- upper.tri(transition_matrix)[idx]
  
  surv_vals <- transition_matrix[idx[!is_recruitment]]
  fec_vals <- transition_matrix[idx[is_recruitment]]
  
  dim <- nrow(transition_matrix)

  fun <- function (landscape, timestep) {
    
    # pull out or create survival/fecundity multipliers
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
    
    
    # get per-cell versions of fecundity and survival values
    survs <- kronecker(surv_mult, t(surv_vals), "*")
    fecs <- kronecker(fec_mult, t(fec_vals), "*")
    
    # empty transition array to fill
    n_cells <- length(cell_idx)
    transition_array <- array(0, c(dim, dim, n_cells))
    
    # convert index from matrix to array
    addition <- dim ^ 2 * (seq_len(n_cells) - 1)
    idx_full <- as.numeric(outer(idx, addition, FUN = "+"))
    
    # put the surv/fec values back in (is_recruitment is recycled to match length)
    transition_array[idx_full[!is_recruitment]] <- survs
    transition_array[idx_full[is_recruitment]] <- fecs
    
    transition_array
    
  }
  
  as.transition_function(fun)
  
}


#' @rdname transition_function
#'
#' @param stages which life-stages contribute to density dependence - default is all
#' @param mask a matrix of boolean values (TRUE/FALSE), equal in dimensions to the life-stage transition matrix
#' and specifying which vital rates (i.e. survival and fecundity) are to be modified by the function
#' @param R_max optional value of maximum growth rate (lambda) if known
#' @param initial_stages optional vector of stable age distributions if known
#' 
#' @export
#'
#' @examples
#' 
#' test_comp_transition <- competition_density(egk_mat)

competition_density <- function(transition_matrix,
                                stages = NULL,
                                mask = NULL,
                                R_max = NULL,
                                initial_stages = NULL) {
  
  dim <- nrow(transition_matrix)

  fun <- function (landscape, timestep) {

    # get metrics and constructor info
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
    n_cells <- length(cell_idx)
    
    # get population matrix
    pop_raster <- landscape$population
    population <- raster::extract(pop_raster, cell_idx)
    
    # get carrying capacity (internal function to STEPS)
    cc <- get_carrying_capacity(landscape, timestep)
    K <- raster::extract(cc, cell_idx)
    
    if (!is.null(stages)) {
      N <- rowSums(population[, stages])
    } else {
      N <- rowSums(population)
    }

    # initialise an array for all of the populations
    transition_array <- array(transition_matrix, dim = c(dim, dim, n_cells))
    
    cells_over <- which(N - K > 0)
    
    # modify life-stage transition matrix and add to array
    for (i in cells_over) {
      transition_array[, , i] <- density_modified_transition(N = N[i],
                                                             K = K[i],
                                                             transition_matrix = transition_matrix,
                                                             mask = mask)
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

get_R <- function (transition_matrix, n_stages = ncol(transition_matrix), initial_stages = NULL, tolerance = 0.001, max_iter = 100) {
  
  if (is.null(initial_stages)) {
    initial_stages <- rep(1, n_stages)
  }
  
  old_stages <- initial_stages
  converged <- FALSE
  iter <- 0
  old_Rs <- rep(0, n_stages)
  
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
  num <- R_max * K
  denom <- R_max * N - N + K
  num / denom
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
find_m <- function(R_target, transition_matrix, mask = NULL, n_stages = ncol(transition_matrix), initial_stages = NULL) {
  
  obj <- function (m, R_target, transition_matrix, mask = NULL, n_stages = ncol(transition_matrix), initial_stages = NULL) {
    new_transition_matrix <- apply_m(m, transition_matrix, mask)
    R_current <- get_R(new_transition_matrix, n_stages = n_stages, initial_stages = initial_stages)
    (R_current - R_target) ^ 2
  } 
  
  out <- stats::optimise(obj, c(0, 1), R_target, transition_matrix, mask, n_stages = n_stages, initial_stages)
  out$minimum
  
}

density_modified_transition <- function (N,
                                         K,
                                         transition_matrix,
                                         n_stages = ncol(transition_matrix),
                                         R_max = NULL,
                                         initial_stages = NULL,
                                         mask = NULL) {
  
  # if the optimal R isn't provided, recalculate it (ideally pre-calculate it to
  # save computation)
  if (is.null(R_max)) {
    R_max <- get_R(transition_matrix, n_stages = n_stages, initial_stages = initial_stages)
  }
  
  # get the target value of R, for this degree of over/under-population
  R_target <- ideal_R(K, N, R_max)
  
  # find a value of m with which to modify transition_matrix, to get to this target value of R
  m <- find_m(R_target, transition_matrix, mask, initial_stages = initial_stages)
  
  # multiply m by the relevant bits of transition_matrix
  transition_matrix_new <- apply_m(m, transition_matrix, mask)
  
  transition_matrix_new
  
}
