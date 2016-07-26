#' @title simulation
#' @rdname simulation
#' @name simulation
#' @description The over-arching function in dlmpr for simulating a-spatial or spatial demographic population projections through time.
#' @param x a dynamic object see \link[dhmpr]{dynamic}
#' @param reps int number of repetitions
#' @param times int number of time steps.
#' @export

simulation <- function(x, reps, times, ...){ 
  
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
  # - are there new patches created? 
  # - are old patches destroyed?
  #
  #  dispersal after manipulation (maybe make it density dependent or based on catastrophy (eg., patch burnt by fire))
  # - re-assess patches 
  # - capture results of patches
  
  #return results
  #summarise simulation (per replication, time-step and overall)
  return(results)
}

#' @rdname simulation
#' @export
is.simulation <- function (x) inherits(x, 'simulation')

# functions to flatten and unflatten population
pop2vec <- function (population) {
  # convert a population dataframe into a vector for deterministic analysis
  ans <- as.vector(t(as.matrix(population)))
  return (ans)
}

vec2pop <- function (vector, population) {
  n_state <- ncol(population)
  n_patch <- nrow(population)
  population[] <- matrix(vector,
                         nrow = n_patch,
                         ncol = n_state,
                         byrow = TRUE)
  return (population)
}

popvecNames <- function (population) {
  # get appropriate names for flattened version of population dataframe
  states <- colnames(population)
  patches <- as.character(seq_len(nrow(population)))
  if (length(patches) == 1) {
    # if only one patch, don't pollute the names
    names <- states
  } else {
    names <- apply(expand.grid(states, patches),
                   1,
                   paste,
                   sep = '',
                   collapse = '_patch_')
  }
  return (names)
}