#' @title habitat_dynamics objects
#' @name habitat_dynamics
#' @rdname habitat_dynamics
#' @description habitat_dynamics are functions for altering the underlying \link[dhmper]{habitat}\code{habitat} at anytime step. habitat_dynamics has to be a function which explicitly changes the habitat. Examples could include: models of fire spread, management actions like culling, the distribution of trawling. 
#' 
#' habitat_dynamics sets up internal or custom functions to work with \code{habitat} and \code{experiment} objects.
#' 
#' @export

as.habitat_dynamics <- function(fun, params, check=FALSE, ...){
  if(!is.function(fun))stop("habitat_dynamics needs to be a function - see the documents for details")
  if(check) {
    message('checking to see your function works with habitat_suitabilit(habitat)')
    test <- do.call(fun,params)
    attr(test, "habitat") <- "habitat"
    return(test)
  }

    fun_params <- structure(list(fun,params),class='habitat_dynamics')
  return(fun_params)
}

#' @rdname habitat_dynamics
#' @export
is.habitat_dynamics <- function (x) inherits(x, 'habitat_dynamics')

#' @rdname habitat_dynamics
#' @name run_habitat_dynamics
#' @export
#' @description this bad boy will run the habitat_dynamics in a experiment.
run_habitat_dynamics <- function(habitat_dynamics, ...){
  if(!is.habitat_dynamics(habitat_dynamics))
  stop("you need to define a habitat_dynamics module in order to run it within an experiment - see the documents for details")
  fun <- habitat_dynamics[[1]]
  params <- habitat_dynamics[[2]]
  altered_habitat <- do.call(fun,params)
  attr(altered_habitat, "habitat") <- "habitat"
  return(altered_habitat)
}  
  