#' How the population changes in a landscape.
#'
#' Pre-defined functions to define population change (e.g. growth) during a simulation.
#'
#' @name population_change_functions
#'
#' @param transition_matrix A symmetrical age-based (Leslie) or stage-based
#'   population structure matrix.
#' @param demographic_stochasticity should demographic stochasticity be used in
#'   population change?
#' @param global_stochasticity,local_stochasticity either scalar values or
#'   matrices (with the same dimension as \code{transition_matrix}) giving
#'   variability (in standard deviations) in the transition matrix either for
#'   all populations (\code{global_stochasticity}) or for each population
#'   separately (\code{local_stochasticity})
#'
#' @rdname population_change_functions
#' 
#' @export
#' 
#' @examples
#' 
#' # Use a function to modify the  
#' # population using life-stage transitions:
#'
#' test_growth <- growth(egk_mat)
growth <- function (transition_matrix,
                    demographic_stochasticity = TRUE,
                    global_stochasticity = 0,
                    local_stochasticity = 0) {
  
  idx <- which(transition_matrix != 0)
  is_recruitment <- upper.tri(transition_matrix)[idx]
  upper <- ifelse(is_recruitment, Inf, 1)
  vals <- transition_matrix[idx]
  dim <- nrow(transition_matrix)
  
  if (is.matrix(global_stochasticity)) {
    stopifnot(identical(dim(transition_matrix), dim(global_stochasticity)))
    stopifnot(identical(which(global_stochasticity != 0), idx))
    global_stochasticity <- global_stochasticity[idx]
  }

  if (is.matrix(local_stochasticity)) {
    stopifnot(identical(dim(transition_matrix), dim(local_stochasticity)))
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
    
    values <- vals + total_noise
    values <- pmax(values, 0)
    values <- pmin(values, rep(upper, n_cells))
    
    transition_array <- array(0, c(dim, dim, n_cells))
    
    # pad the index to get corresponding elements in each slice  
    addition <- length(transition_matrix) * (seq_len(n_cells) - 1)
    idx_full <- as.numeric(outer(idx, addition, FUN = "+"))
    transition_array[idx_full] <- values

    if (demographic_stochasticity) {
      
      local_t <- array(apply(transition_array, 3,
                             function(x) {
                               rbind(0,
                                     x[-1,],
                                     if (is.vector(x[-1, ])) {
                                       x[-1,]
                                     } else {
                                       1 - apply(x[-1, ], 2, sum)
                                     })
                             }),
                       dim = c(dim + 1,
                               dim,
                               n_cells))
      
      local_f <- transition_array
      local_f[-1, , ] <- 0
      
      pop_tmp <- cbind(population, rep(0, n_cells))
      
      survival_stochastic <- sapply(seq_len(ncol(population)),
                                    function(y) sapply(seq_len(nrow(population)),
                                                       function(x) stats::rmultinom(n = 1,
                                                                                    size = pop_tmp[x, y],
                                                                                    prob = local_t[, y, x])),
                                    simplify = 'array')
      
      
      new_offspring_deterministic <- sapply(seq_len(nrow(population)), function(x) local_f[ , , x] %*% matrix(population[x, ]))
      new_offspring_stochastic <- matrix(stats::rpois(n = length(c(new_offspring_deterministic)),
                                                      lambda = c(new_offspring_deterministic)),
                                         nrow = nrow(new_offspring_deterministic))
      new_offspring <- apply(new_offspring_stochastic, 2, sum)
      
      population <- t(apply(survival_stochastic[seq_len(ncol(population)), , ], c(1, 2), sum))
      population[ , 1] <- population[ , 1] + new_offspring
      
    } else {
      
      population <- t(sapply(seq_len(n_cells),
                             function(x) transition_array[ , , x] %*% matrix(population[x, ])))
      
      # get whole integers
      population <- t(apply(population,
                            1,
                            FUN = function(x) rbinom(prob = (x/ceiling(max(population))),
                                                     size = ceiling(max(population)),
                                                     n = 2)))
    }
    
    
    # put back in the raster
    population_raster[cell_idx] <- population
    
    landscape$population <- population_raster
    
    landscape
  }
  
  as.population_growth(pop_dynamics)
  
}


# spatial_transition <- function (survival_raster,
#                                 fecundity_raster,
#                                 demographic_stochasticity = TRUE,
#                                 global_stochasticity = 0,
#                                 local_stochasticity = 0) {
#   
#   # as above
#   
# }


# functional_transition <- function (survival_function,
#                                    fecundity_function,
#                                    demographic_stochasticity = TRUE,
#                                    global_stochasticity = 0,
#                                    local_stochasticity = 0) {
#   
#   # as above
#   
# }

# # e.g.:
# survival_fun <- function (landscape) {
#   cov <- landscape$time_since_fire
#   if (is.null(suit)) stop ()
#   
#   cell_idx <- which(!is.na(raster::getValues(cov)))
#   cov_vec <- raster::extract(cov, cell_idx)
#   
#   multipliers <- c(0.5, 0.3, 0.2)
#   cov_mat <- outer(log(cov_vec), multipliers, FUN = "*")
#   plogis(cov_mat)
# 
# }


# fecundity_fun <- function (landscape) {
#   cov <- landscape$time_since_fire
#   if (is.null(suit)) stop ()
#   
#   cell_idx <- which(!is.na(raster::getValues(cov)))
#   cov_vec <- raster::extract(cov, cell_idx)
#   
#   fec_vec <- exp(log(cov_vec) * 0.8)
#   cbind(0, 0, fec_vec)
# }

##########################
### internal functions ###
##########################

as.population_growth <- function (simple_growth) {
  as_class(simple_growth, "population_dynamics", "function")
}