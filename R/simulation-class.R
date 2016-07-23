#' @title simulation
#' @rdname simulation
#' @name simulation
#' @description The over-arching function in dlmpr for simulating a-spatial or spatial demographic population projections through time.
#' @param dots that can contain:
#'  transition a transition object, See \link[dlmpr]{as.transition}.
#'  population a population object, see \link[dlmpr]{as.population}
#'  habitat a habitat object, see \link[dlmpr]{as.habitat}.
#'  dispersal a dispersal object, see \link[dlmpr]{as.dispersal}. 
#'  module a module object, see \link[dlmpr]{as.module}.
#' @export

simulation <- function(...){ #transition,population,habitat,dispersal,module){
  #introduce transtion, do checks.
  
  #introduce population, do checks.
  
  #introduce habitat, do checks.
  
  #introduce dispersal, do checks.
  
  #introduce module(s), do checks.
  
  # capture dhmpr objects
  object <- captureDots(...)
  
  #demographic growth across a-spatial or spatial population(s).
  
  #population management and/or habitat alteration (via module(s)).
  #need to develop a function that updates habitat
  
  #dispersal(s) or stages or stage after module (if any).
  
  #population management after dispersal.
  
  #summarise simulation (per replication, time-step and overall)
  return(object)
}

#' @rdname simulation
#' @export
is.simulation <- function (x) inherits(x, 'simulation')

captureDots <- function (...) {
  # capture arguments passed as dots, and grab their names even if not directly
  # named
  ans <- list(...)
  dots <- substitute(list(...))[-1]
  names(ans) <- sapply(dots, deparse)
  ans
}

