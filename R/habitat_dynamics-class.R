#' @title habitat_suitability_dynamics objects
#' @name habitat_suitability_dynamics
#' @rdname habitat_suitability_dynamics
#' @description habitat_suitability_dynamics are functions for altering the underlying \link[dhmper]{habitat}\code{habitat_suitability} at anytime step. Habitat_suitability_dynamics has to be a function which explicitly changes the habitat_suitability. Examples could include: models of fire spread, management actions like culling, the distribution of trawling. 
#' 
#' habitat_suitability_dynamics sets up internal or custom functions to work with \code{habitat} and \code{experiment} objects.
#' 
#' @export
#' @examples 
#' ## Create population
#' library(raster)
#' library(dhmpr)
#' #set a transition matrix
#' mat <- matrix(c(.53,0,.62,0.15,0.87,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult') 
#' trans <- as.transition(mat)
#' n_stages <- length(states(trans))
#' 
#' set up starting populations
#' set.seed(42)
#' xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
#' Dd <- as.matrix(dist(xy))
#' w <- exp(-1/nrow(xy) * Dd)
#' Ww <- chol(w)
#' xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
#' coordinates(xy) <- ~x+y
#' r <- rasterize(xy, raster(points2grid(xy)), 'z')
#' proj4string(r) <- '+init=epsg:4283'
#' r[] <- scales::rescale(r[],to=c(0,1))
#' random_populations <- sampleRandom(r, size=50, na.rm=TRUE, sp=TRUE)
#' random_populations@data <- as.data.frame(t(rmultinom(50, size = 50, prob = c(0.8,0.2,0.1))))
#' 
#' # set up habitat
#' features <- list('habitat_suitability_map'=as.habitat_suitability(r),
#'                  'population'=as.populations(random_populations),
#'                  'carrying_capacity'=as.carrying_capacity(100))
#' habitat <- as.habitat(features)
#' 
#' ##Define a function for manipulating habitat.
#' ##This can be a custom function for manipulating rasters or existing functions. 
#' ##Create a named list with corresponding parameters and values
#' params = list(habitat,
#'              fire_start_location = sample(ncell(habitat_suitability(habitat)),10),
#'              prob = 0.24,
#'              continue_to_burn_prob = 0.01)
#'               
#'## check habitat_suitability_dynamics produces expected output
#' fire_affected_habitat_suitability <- as.habitat_suitability_dynamics(fire_module,params,check=TRUE)                   
#'
#'## If it does? create the habitat_suitability_dynamics (yay!).
#'fire_habitat_suitability_dynamics <- as.habitat_suitability_dynamics(fun,params) 
#'
as.habitat_suitability_dynamics <- function(fun, params, check=FALSE, ...){
  if(!is.function(fun))stop("habitat_suitability_dynamics needs to be a function - see the documents for details")
  if(check) {
    message('checking to see your function works with habitat_suitabilit(habitat)')
    test <- do.call(fun,params)
    attr(test, "habitat") <- "habitat_suitability"
    return(test)
  }

    fun_params <- structure(list(fun,params),class='habitat_suitability_dynamics')
  return(fun_params)
}

#' @rdname habitat_suitability_dynamics
#' @export
is.habitat_suitability_dynamics <- function (x) inherits(x, 'habitat_suitability_dynamics')
#' @rdname run_habitat_suitability_dynamics
#' @export
#' @description this bad boy will run the habitat_suitability_dynamics in a experiment.
run_habitat_suitability_dynamics <- function(habitat_suitability_dynamics, ...){
  if(!is.habitat_suitability_dynamics(habitat_suitability_dynamics))
  stop("you need to define a habitat_suitability_dynamics module in order to run it within an experiment - see the documents for details")
  fun <- habitat_suitability_dynamics[[1]]
  params <- habitat_suitability_dynamics[[2]]
  altered_habitat_suitability <- do.call(fun,params)
  attr(altered_habitat_suitability, "habitat") <- "habitat_suitability"
  return(altered_habitat_suitability)
}  
  