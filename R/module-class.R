#' @title module objects
#' @name module
#' @rdname module
#' @description Module for altering the habitat, a model could be fire spread, management act, the distribution of trawling. Module sets up internal or custom functions to work with \code{habitat} and \code{experiment} objects.
#' 
#' Module needs to be a function, that gets imported into the simulaion function. 

#' @export
#' @examples 
#' ## Create population
#' pop <- as.population(data.frame('larvae'=80,'juvenile'=29,'adult'=5) )
#'
#'##Set up habitat.
#' library(raster)
#' set.seed(42)
#' xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
#' Dd <- as.matrix(dist(xy))
#' w <- exp(-1/nrow(xy) * Dd)
#' Ww <- chol(w)
#' xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
#' coordinates(xy) <- ~x+y
#' r <- rasterize(xy, raster(points2grid(xy)), 'z')
#' hab <- raster(r)
#' res(hab) <- 0.01
#' hab <- resample(r, hab)
#' proj4string(hab) <- '+init=epsg:4283'
#' habs <- as.habitat(list(hab,population = pop))
#'
#'##Define a function for manipulating habitat.
#'##This can be a custom function for manipulating rasters or existing functions. 
#'fun <- fire_spread
#'
#'##Create a named list with corresponding parameters and values
#'params = list(habitat=habs,
#'              fire_start_location = sample(ncell(suitability(habs)),10),
#'              prob = 0.24,
#'              continue_to_burn_prob = 0.01)
#'               
#'## check module produces expected output
#'fire_module <- as.module(fun,params,check=TRUE)                             
#'
#'## If it does? create the module (yay!).
#'fire_module <- as.module(fun,params) 
#'

as.module <- function(fun, params, check=FALSE, ...){
  
  stopifnot(is.function(fun))
  if(check) {
    test <- do.call(fun,params)
    return(test)
  }
  
  args <- names(formals(fun))
  if (args[1] != 'habitat') {
    stop ("module must contain 'habitat' as the first parameter")
  }
  
  attr(fun, 'user-defined-module') <- TRUE
  fun_params <- structure(list(fun,params),class='module')

    return(fun_params)

}

#' @rdname module
#' @export
is.module <- function (x) inherits(x, 'module')
