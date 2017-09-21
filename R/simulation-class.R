#' @title simulation
#' @rdname simulation
#' @name simulation
#' @description The over-arching function in dhmpr for simulating cell-based spatial demographic meta-population projections through time in a dynamic landscape (or seascape).
#' @param x an dynamic habitat metapopulation experiment object see \link[dhmpr]{experiment}
#' @param reps int number of simulations to run.
#' @param times int The number of time steps of the least frequent discrete event. e.g. yearly population growth/management.
#' @export

simulation <- function(x, reps, ...){ 
  
  # I have created a dynamics object which will contain all the relevant info
  # capture dhmpr objects
  hab <- habitat(x)
  pops <- pop2vec(x[['population']])
  proj_mat <- as.matrix(x)
  
  if(nrow(coords(hab))==1){
  results <- demographic(as.population(pops),as.transition(proj_mat),nrep = reps,
              time = times)
  class(results) <- "demographic"
  return(results)
  }  
  # ideally we can do this with C++ - to make it fast.
  #
  # - We'll need an input of A = the projection matrix
  # - population
  # - keep record of patches
  # - module manipualted habitat 
  # - are there new patches created? _ implenet fire model across habitat/ land scape distrubance/ management.  
  # - are old patches destroyed?
  #
  #  dispersal after manipulation (maybe make it density dependent or based on catastrophy (eg., patch burnt by fire))
  # - run FFT here. Recalulate populations. 
  # - re-assess patches 
  # - capture results of patches
  
  #return results
  #summarise simulation (per replication, time-step and overall)
  return(results)
}

#' @rdname simulation
#' @export
is.simulation <- function (x) inherits(x, 'simulation')