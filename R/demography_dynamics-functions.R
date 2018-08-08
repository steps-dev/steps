#' Functions to modify the demography in a state object.
#' 
#' Pre-defined functions to operate on a population demography
#' during a simulation.
#'
#' @name demography_dynamics_functions
#'
#' @param transition_matrix A life-stage transition matrix.
#' @param stochasticity A matrix with standard deviations (consistent or
#' varying) around the transition means with dimensions matched to the
#' life-stage transition matrix or a number representing a consistent
#' standard deviation to apply to all transitions (default is 0).
#' @param fecundity_fraction A multiplier value between 0 and 1 for fecundity
#' values in the transition matrix.
#' @param survival_fraction A multiplier value between 0 and 1 for survival
#' values in the transition matrix.
#' @param surv_layers a list of raster stacks with multipliers for survival
#' equal to the number of life-stages.
#' @param fec_layers a list of raster stacks with multipliers for fecundities
#' equal to the number of life-stages. Note, life-stages that do not reproduce
#' will have NULL in place of the raster stack.
#' 
#' @examples
#' 
#' library(steps)
#' 
#' mat <- matrix(c(0.000,0.000,0.302,0.302,
#'                 0.940,0.000,0.000,0.000,
#'                 0.000,0.884,0.000,0.000,
#'                 0.000,0.000,0.793,0.793),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#' 
#' mat_sd <- .01

#' @rdname demography_dynamics_functions
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the demo_environmental_stochasticity function to modify the transition
#' # matrix with specified environmental stochasticity:
#' test_demo_es <- demo_environmental_stochasticity(transition_matrix = mat,
#'                                     stochasticity = mat_sd)

demo_environmental_stochasticity <- function (transition_matrix,
                                              stochasticity=0) {
  
  idx <- which(transition_matrix != 0)
  is_recruitment <- upper.tri(transition_matrix)[idx]
  lower <- 0
  upper <- ifelse(is_recruitment, Inf, 1)
  vals <- transition_matrix[idx]
  
  if (is.matrix(stochasticity)) {
    stopifnot(identical(dim(transition_matrix), dim(stochasticity)))
    stopifnot(identical(which(stochasticity != 0), idx))
    stochasticity <- stochasticity[idx]
  }
  
  env_stoch_fun <- function (state, timestep) {
    
    # if (!is.null(state$demography$local_transition_matrix) & is.null(local_transition_matrix)) {
    #   stop("Local cell-based transition matrices are required \nfor this function - none have been specified")
    # }
    # global_demography <- transition_matrix
    # nstages <- dim(transition_matrix)[[1]]
    
    # demography_obj <- state$demography$demography_obj
    
    local <- !is.null(state$demography$local_transition_matrix)
    
    if (local) {
      demography_obj <- state$demography$local_transition_matrix
    } else {
      demography_obj <- state$demography$global_transition_matrix
    }
    
    # change to:
    # demography_obj <- state$demography$demography_obj
    
    idx <- replicate_values(idx, demography_obj, index = TRUE)
    vals <- replicate_values(vals, demography_obj)
    lower <- replicate_values(lower, demography_obj)
    upper <- replicate_values(upper, demography_obj)
    stochasticity <- replicate_values(stochasticity, demography_obj)
    
    demography_obj[idx] <- extraDistr::rtnorm(length(idx),
                                              vals,
                                              stochasticity,
                                              a = lower,
                                              b = upper)
    
    if (local) {
      state$demography$local_transition_matrix <- demography_obj
    } else {
      state$demography$global_transition_matrix <- demography_obj
    }
    
    # change to:
    # state$demography$demography_obj <- demography_obj
    
    state
    
  }
  
  as.demography_environmental_stochasticity(env_stoch_fun)
  
}


#' @rdname demography_dynamics_functions
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the demo_density_dependence function to modify the transition
#' # matrix once carrying capacity is reached:
#' test_demo_dd <- demo_density_dependence(transition_matrix = mat,
#'                                         fecundity_fraction = 1,
#'                                         survival_fraction = 0.5)

demo_density_dependence <- function (transition_matrix,
                                     fecundity_fraction = 1,
                                     survival_fraction = 1) {
  
  idm <- which(transition_matrix != 0)
  fecundity <- which(transition_matrix != 0 & upper.tri(transition_matrix))
  survival <- setdiff(idm, fecundity)
  vals_fecundity <- transition_matrix[fecundity]
  vals_survival <- transition_matrix[survival]
  
  dens_dep_fun <- function (state, timestep) {
    
    idr <- which(!is.na(raster::getValues(state$population$population_raster[[1]])))
    population <- raster::extract(state$population$population_raster, idr)
    
    local <- !is.null(state$demography$local_transition_matrix)
    
    if (local) {
      demography_obj <- state$demography$local_transition_matrix
    } else {
      demography_obj <- state$demography$global_transition_matrix
    }
    
    if (any(rowSums(population) > raster::extract(state$habitat$carrying_capacity, idr))) {
      
      # identify cells that are still within carrying capacity
      idk <- which(rowSums(population) <= raster::extract(state$habitat$carrying_capacity, idr))
      
      # modify transition matrices
      vals_fecundity <- replicate_values(vals_fecundity, demography_obj)
      vals_survival <- replicate_values(vals_survival, demography_obj)
      fecundity <- replicate_values(fecundity, demography_obj, index = TRUE)
      survival <- replicate_values(survival, demography_obj, index = TRUE)
      
      demography_obj[fecundity] <- vals_fecundity * fecundity_fraction
      demography_obj[survival] <- vals_survival * survival_fraction
      
      # if local transition matrices are used, restore original values for
      # cells where population is within carrying capacity
      if (local) {
        demography_obj[, , idk] <- transition_matrix
      }
      
    }
    
    if (local) {
      state$demography$local_transition_matrix <- demography_obj
    } else {
      state$demography$global_transition_matrix <- demography_obj
    }
    
    state
    
  }
  
  as.demography_density_dependence(dens_dep_fun)
  
}


#' @rdname demography_dynamics_functions
#'
#' @export
#' 
#' @examples
#' 
#' # Use the demo_surv_fec_modify function to modify the  
#' # demography using explicit survival and fecundity layers:
#' test_survfec <- demo_surv_fec_modify(transition_matrix = mat,
#'                                     surv_layers = surv,
#'                                     fec_layers = fec)

demo_surv_fec_modify <- function (transition_matrix, surv_layers, fec_layers) {
  
  surv_fec_mod <- function (state, timestep) {
    
    if (is.null(state$demography$local_transition_matrix)) {
      stop("Local cell-based transition matrices are required \nfor this function - none have been specified")
    }
    
    global_demography <- transition_matrix
    nstages <- dim(transition_matrix)[[1]]
    
    if (any(lapply(surv_layers, function (x) is.null(x)) == TRUE)) {
      stop("Survival layers must be provided for all life-stages")
    }
    
    if (any(unlist(lapply(surv_layers, function (x) raster::nlayers(x))) < timestep)) {
      stop("The number of survival/fecundity layers must match \nthe number of timesteps in the simulation run")
    }
    
    local_demography <- state$demography$local_transition_matrix
    
    for (i in seq_len(nstages)) {
      
      matrix_idx <- which(global_demography != 0 & row(global_demography) != 1 & col(global_demography) == i, arr.ind = TRUE)
      local_demography[matrix_idx[1], matrix_idx[2], ] <- global_demography[matrix_idx] * surv_layers[[i]][[timestep]][]
    }
    
    for (i in seq_len(nstages)) {
      
      if (!is.null(fec_layers[[i]])) {
        local_demography[1, i, ] <- global_demography[1, i] * fec_layers[[i]][[timestep]][]
      }
      
    }
    
    state$demography$local_transition_matrix <- local_demography
    
    state
    
  }
  
  as.demography_modify_surv_fec(surv_fec_mod)
  
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

as.demography_modify_surv_fec <- function (demography_modify_surv_fec) {
  as_class(demography_modify_surv_fec, "demography_dynamics", "function")
}

