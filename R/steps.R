#' @title Simulate population trajectories over space and time with custom
#' dynamic functions.
#' @name steps
#' @docType package
#' @description Simulating shifts in species populations is an important
#' part of ecological management. Species respond to spatial and temporal
#' changes in the landscape resulting from environmental phenomena, managerial
#' actions or anthropogenic activities. This data is crucial for modelling,
#' however, current software that incorporates this information has limited
#' flexibility, transparency, and availability. \code{steps} extends the features
#' found in existing software and uses common spatial inputs that are derived
#' from many other software packages.
#'
#' A \link[steps]{simulation} is run on a \link[steps]{landscape} using population
#' dynamics functions contained in a \link[steps]{population_dynamics} object.
#' \link[steps]{habitat_dynamics_functions} can also be added to the simulation to
#' modify the habitat during a simulation.

steps_stash <- new.env()

flush_stash <- function() {
  for (name in names(steps_stash)) {
    steps_stash[[name]] <- NULL
  }
}