#' @title module objects
#' @name module
#' @rdname module
#' @description Module for altering the habitat, they maybe fire spread, management, ect. 
#' Module sets up internal or custom functions to work with \code{habitat} and \code{simulation} objects.
#' 
#' Module needs to be a function, that gets imported into the simulaion function. 

#' @export
#' @examples 
#'
#' params = list(habitat=habitat,
#'               fire_start_location = sample(ncell(r2),10),
#'               prob = 0.24,
#'               continue_to_burn_prob = 0.01,
#'               id=TRUE)
#'               
#' fun <- fire_spread                


as.module <- function(fun, params, ...){
  
  stopifnot(is.function(fun))
  params <- params
  environment(fun) <- environment()
  args <- names(formals(fun))
  if (args[1] != 'habitat') {
    stop ("module must contain 'habitat' as the first parameter")
  }
  
  attr(fun, 'user-defined-module') <- TRUE
  
  # test it runs.
  # do.call(fun,params) 
  
  return(list(fun,params))
  
}
