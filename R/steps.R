#' @title Simulate population trajectories over space and time with dynamic functions.
#' @name steps
#' @aliases steps-package
#' @description Simulating shifts in species populations is an important
#' part of ecological management. Species respond to spatial and temporal
#' changes in the landscape resulting from environmental phenomena, managerial
#' actions or anthropogenic activities. This data is crucial for modelling,
#' however, current software that incorporates this information has limited
#' flexibility, transparency, and availability. \code{steps} extends the features
#' found in existing software and accepts common spatial inputs that are derived
#' from many other existing software packages.
#'
#' A \link[steps]{simulation} is run on a \link[steps]{landscape} using population
#' dynamics functions contained in a \link[steps]{population_dynamics} object.
#' \link[steps]{habitat_dynamics_functions} can also be added to the simulation to
#' modify the habitat during a simulation.
"_PACKAGE"


if (!exists("steps_stash", envir = .GlobalEnv)) {
  steps_stash <- new.env()
}
# steps_stash environment- only create if does not exist
# KH - no longer supporting multisessions because of terra pointer issue
# without this change, simulations crash because environments continue to renew and terra pointers are lost

flush_stash <- function() {
  for (name in names(steps_stash)) {
    steps_stash[[name]] <- NULL
  }
}

# replace the values in the steps stash with those in this new one (used to pass the stash onto parallel workers)
replace_stash <- function(new_stash) {

  # flush the old one
  flush_stash()
  
  for (name in names(new_stash)) {
    steps_stash[[name]] <- new_stash[[name]]
  }
  
}

options(future.globals.maxSize= 1000*1024^2)
terra::terraOptions(memmax = 1e+09)
