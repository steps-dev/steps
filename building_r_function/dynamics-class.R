#' @title dynamics objects
#' @name dynamics
#' @rdname dynamics
#' @description dynamics are functions for altering the underlying \link[dhmper]{habitat}\code{habitat} at anytime step. dynamics has to be a function which explicitly changes the habitat. Examples could include: models of fire spread, management actions like culling, the distribution of trawling. 
#' 
#' dynamics sets up internal or custom functions to work with \code{habitat} and \code{experiment} objects.
#' 
#' @export

as.dynamics <- function(fun, params, check=FALSE, ...){
  if(!is.function(fun))stop("dynamics needs to be a function - see the documents for details")
  if(check) {
    message('checking to see your function works with habitat_suitabilit(habitat)')
    test <- do.call(fun,params)
    attr(test, "habitat") <- "habitat"
    return(test)
  }

    fun_params <- structure(list(fun,params),class='dynamics')
  return(fun_params)
}

#' @rdname dynamics
#' @export
is.dynamics <- function (x) inherits(x, 'dynamics')

#' @rdname dynamics
#' @name run_dynamics
#' @export
#' @description this bad boy will run the dynamics in a experiment.
run_dynamics <- function(dynamics, habitat_object, time_step){
  if(!is.dynamics(dynamics))
  stop("you need to define a dynamics module in order to run it within an experiment - see the documents for details")
  fun <- dynamics[[1]]
  
  if(inherits(dynamics[[2]],c("RasterLayer"))){
    params <- list(habitat_object, dynamics[[2]])
  } else {
    params <- list(habitat_object, dynamics[[2]][[time_step]])
  }
  
  altered_habitat <- do.call(fun,params)
  attr(altered_habitat, "habitat") <- "habitat"
  return(altered_habitat)
}  

# #' @rdname dynamics
# #' @name run_dynamics
# #' @export
# #' @description this bad boy will run the dynamics in a experiment.
# run_dynamics <- function(dynamics, ...){
#   if(!is.dynamics(dynamics))
#     stop("you need to define a dynamics module in order to run it within an experiment - see the documents for details")
#   fun <- dynamics[[1]]
#   params <- dynamics[[2]]
#   altered_habitat <- do.call(fun,params)
#   attr(altered_habitat, "habitat") <- "habitat"
#   return(altered_habitat)
# } 
