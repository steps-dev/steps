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
# @param env_stoch a function for adding environmental stochasticity to the global transition matrix at each timestep
# @param determ_surv_fec a function for altering a life-stage transition matrix with user supplied spatial layers at each timestep
# @param demo_dens_dep a function for modifying the transition matrix at each timestep when carrying capacity is reached
#' @param transition_matrix a life-stage transition matrix
#' @param stochasticity a matrix with standard deviations (consistent or varying) around the transition means with matching dimensions as the life-stage transition matrix or a number representing a consitent standard deviation to apply to all transitions (default is 0)
#' @param fecundity_fraction a multiplier value between 0 and 1 for fecundity values in the transition matrix
#' @param survival_fraction a multiplier value between 0 and 1 for fecundity values in the transition matrix 
#' @param surv_layers a list of raster stacks with multipliers for survival equal to the number of life-stages
#' @param fec_layers a list of raster stacks with multipliers for fecundities equal to the number of life-stages. Note, life-stages that do not reproduce will have NULL in place of the raster stack
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
#' r2 <- r
#' r2[na.omit(r2)] <- sample(r[na.omit(r)])
#' 
#' r3 <- r
#' r3[na.omit(r3)] <- sample(r[na.omit(r)])
#' 
#' surv <- list(stack(r2, r2, r2),
#'              stack(r2, r2, r2),
#'              stack(r2, r2, r2),
#'              stack(r2, r2, r2))
#'              
#' fec <- list(NULL,
#'             NULL,
#'             stack(r3, r3, r3),
#'             stack(r3, r3, r3))
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
#'   state
#' }
#' 
#' # Define the function as a demography_dynamics object
#' 
#' example_function <- as.demography_dynamics(example_function)

as.demography_dynamics <- function (demography_dynamics_function) {
  as_class(demography_dynamics_function, "demography_dynamics", "function")
}

#' @rdname demography_dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'demography dynamics'
#'   
#' is.demography_dynamics(example_function)

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
#' print(example_function)

print.demography_dynamics <- function (x, ...) {
  cat("This is a demography_dynamics object")
}


#' @rdname demography_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the demography_dynamics function to modify a demography object:
#' 
#' env_stoch <- demo_environmental_stochasticity(transition_matrix = mat,
#'                                              stochasticity = mat_sd)
#' 
#' example_function <- demography_dynamics(env_stoch,
#'                                         demo_dens_dep =  demo_density_dependence())

demography_dynamics <- function (...) {
  
  dots <- list(...)
  
  # run checks on the functions they've passed in, make sure they are legit
  
  demo_dynamics <- function (state, timestep) {
    
    if (!is.null(unlist(dots))){
      for (fun in dots) {
        state <- fun(state, timestep)
      }
    }
    
    state
    
  }
  
  as.demography_dynamics(demo_dynamics)
  
}

# demography_dynamics <- function (env_stoch = NULL,
#                                  determ_surv_fec = NULL,
#                                  demo_dens_dep = NULL) {
#   
#   demographic_dynamics <- function (state, timestep) {
#     
#     if (!is.null(env_stoch))
#       state <- env_stoch(state, timestep)
#     
#     if (!is.null(determ_surv_fec))
#       state <- determ_surv_fec(state, timestep)
#     
#     if (!is.null(demo_dens_dep))
#       state <- demo_dens_dep(state, timestep)
# 
#     state
#     
#   }
#   
#   as.demography_dynamics(demographic_dynamics)
#   
# }

##########################
### internal functions ###
##########################

as.demography_environmental_stochasticity <- function (demography_environmental_stochasticity) {
  as_class(demography_environmental_stochasticity, "demography_dynamics", "function")
}

as.demography_density_dependence <- function (demography_density_dependence) {
  as_class(demography_density_dependence, "demography_dynamics", "function")
}

as.demography_deterministic_surv_fec <- function (demography_deterministic_surv_fec) {
  as_class(demography_deterministic_surv_fec, "demography_dynamics", "function")
}

####################################
### pre-defined module functions ###
####################################

#' @rdname demography_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the demo_environmental_stochasticity function to modify the transition
#' # matrix with specified environmental stochasticity:
#' 
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
    
    # expand these values if demography_obj is an array (otherwise leave them as they are)
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


#' @rdname demography_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the demo_density_dependence function to modify the transition
#' # matrix once carrying capacity is reached:
#' 
#' test_demo_dd <- demo_density_dependence(fecundity_fraction = 1,
#'                                                 survival_fraction = 0.5)

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

      idk <- which(rowSums(population) <= raster::extract(state$habitat$carrying_capacity, idr))
      
      vals_fecundity <- replicate_values(vals_fecundity, demography_obj)
      vals_survival <- replicate_values(vals_survival, demography_obj)
      fecundity <- replicate_values(fecundity, demography_obj, index = TRUE)
      survival <- replicate_values(survival, demography_obj, index = TRUE)

      demography_obj[fecundity] <- vals_fecundity * fecundity_fraction
      demography_obj[survival] <- vals_survival * survival_fraction
      
      demography_obj[, , idk] <- transition_matrix
      
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


#' @rdname demography_dynamics
#'
#' @export
#' 
#' @examples
#' 
#' # Use the deterministic_surv_fec function to modify the  
#' # demography using explicit survival and fecundity layers:
#' 
#' test_survfec <- deterministic_surv_fec(transition_matrix = mat,
#'                                     surv_layers = surv,
#'                                     fec_layers = fec)

surv_fec_modify <- function (transition_matrix, surv_layers, fec_layers) {
  
  surv_fec_mod <- function (state, timestep) {
    
    if (is.null(state$demography$local_transition_matrix)) {
      stop("Local cell-based transition matrices are required \nfor this function - none have been specified")
    }
    
    global_demography <- transition_matrix
    nstages <- dim(transition_matrix)[[1]]
    
    if (any(unlist(lapply(surv_layers, function (x) raster::nlayers(x))) < timestep)) {
      stop("The number of survival/fecundity layers must match \nthe number of timesteps in the simulation run")
    }
    
    local_demography <- state$demography$local_transition_matrix
    
    for (i in seq_len(nstages)) {
      
      if (any(lapply(surv_layers, function (x) is.null(x)) == TRUE)) {
        stop("Survival layers must be provided for all life-stages")
      }
      
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
  
  as.demography_deterministic_surv_fec(surv_fec_mod)
  
}