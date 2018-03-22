#' @useDynLib steps
#' @importFrom Rcpp sourceCpp
NULL
#' Change the population in a state object
#' 
#' @description A 'population dynamics' object is used to modify species populations in space and time.
#' It is a sub-component of a \link[steps]{dynamics} object and is executed in each timestep of an experiment.
#'
#' @rdname population_dynamics
#'
#' @param population_dynamics_function A function that operates on a state object to change population at specified timesteps. User may enter a custom function or select a pre-defined module - see documentation. 
#' @param x a population_dynamic object
#' @param ... further arguments passed to or from other methods
#' @param state a state object to apply the population function to
#' @param timestep the timestep in the experiment to apply the population function to the state object
#'
#' @return An object of class \code{population_dynamics}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' mat <- matrix(c(0.000,0.000,0.302,0.302,
#'                 0.940,0.000,0.000,0.000,
#'                 0.000,0.884,0.000,0.000,
#'                 0.000,0.000,0.793,0.793),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#' 
#' pop <- stack(replicate(4, ceiling(r * 0.2)))
#' 
#' hab_suit <- r / cellStats(r, "max")
#' 
#' k <- ceiling(hab_suit * 10)
#' 
#' dp_params <- list(dispersal_distance=list('Stage_1'=0,
#'                                           'Stage_2'=10,
#'                                           'Stage_3'=10,
#'                                           'Stage_4'=0),
#'                   dispersal_kernel=list('Stage_1'=0,
#'                                         'Stage_2'=exp(-c(0:9)^1/3.36),
#'                                         'Stage_3'=exp(-c(0:9)^1/3.36),
#'                                         'Stage_4'=0),
#'                   dispersal_proportion=list('Stage_1'=0,
#'                                             'Stage_2'=0.35,
#'                                             'Stage_3'=0.35*0.714,
#'                                             'Stage_4'=0)
#'                   )
#' 
#' test_habitat <- build_habitat(habitat_suitability = hab_suit,
#'                               carrying_capacity = k)
#' test_demography <- build_demography(transition_matrix = mat,
#'                                     dispersal_parameters = rlnorm(1))
#' 
#' test_demography_dp <- build_demography(transition_matrix = mat,
#'                                     dispersal_parameters = dp_params)
#' test_population <- build_population(pop)
#' 
#' test_state <- build_state(test_habitat, test_demography, test_population)
#' 
#' test_state_dp <- build_state(test_habitat,
#'                              test_demography_dp,
#'                              test_population)
#'                              
#' example_function <- function (state, timestep) {
#'  state
#' }
#' 
#' no_population_dynamics <- as.population_dynamics(example_function)

as.population_dynamics <- function (population_dynamics_function) {
  stopifnot(inherits(population_dynamics_function,"function"))
  set_class(population_dynamics_function, "population_dynamics")
}

#' @rdname population_dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'population dynamics'
#'   
#' is.population_dynamics(no_population_dynamics)

is.population_dynamics <- function (x) {
  inherits(x, 'population_dynamics')
}

#' @rdname population_dynamics
#'
#' @export
#'
#' @examples
#' 
#' print(no_population_dynamics)

print.population_dynamics <- function (x, ...) {
  cat("This is a population_dynamics object")
}

##########################
### internal functions ###
##########################

cap_population <- function (new_population, carrying_capacity) {
  
  #carrying_capacity <- state$habitat$carrying_capacity
  
  if (is.null(carrying_capacity)) {
    stop ("carrying capacity must be specified",
          call. = FALSE)
  }
  
  # get degree of overpopulation, and shrink accordingly
  overpopulation <- as.vector(carrying_capacity) / rowSums(new_population)
  overpopulation[is.nan(overpopulation)] <- 0
  overpopulation <- pmin(overpopulation, 1)
  new_population <- sweep(new_population, 1, overpopulation, "*")
  
  new_population
}


dispersal_matrix <- function (locations, distance_decay = 0.5) {
  D <- as.matrix(stats::dist(locations))
  dispersal_matrix <- exp(-D / distance_decay)
  sums <- colSums(dispersal_matrix)
  dispersal_matrix <- sweep(dispersal_matrix, 2, sums, "/")
  dispersal_matrix
}


extend <- function (x, factor = 2) {
  # given an evenly-spaced vector `x` of cell centre locations, extend it to the
  # shortest possible vector that is at least `factor` times longer, has a 
  # length that is a power of 2 and nests the vector `x` approximately in the 
  # middle. This is used to define each dimension of a grid mapped on a torus
  # for which the original vector x is approximately a plane.
  
  # get cell number and width
  n <- length(x)
  width <- x[2] - x[1]
  
  # the smallest integer greater than or equal to than n * factor and an
  # integer power of 2
  n2 <- 2 ^ round(log2(factor * n))
  
  # find how much to pad each end of n
  pad <- n2 - n
  if (pad %% 2 == 0) {
    # if it's even then that's just tickety boo
    pad1 <- pad2 <- pad / 2
  } else {
    # otherwise put the spare cell on the right hand side
    pad1 <- (pad - 1) / 2
    pad2 <- (pad + 1) / 2
  }
  
  # define the padding vectors
  padx1 <- x[1] - rev(seq_len(pad1)) * width 
  padx2 <- x[n] + seq_len(pad2) * width
  
  # combine these and add an attribute returning the start and end indices for
  # the true vector
  newx <- c(padx1, x, padx2)
  attr(newx, 'idx') <- pad1 + c(1, n)
  newx
}


bcb <- function (x, y, f = I) {
  # get a the basis vector for a block-circulant matrix representing the 
  # dispersal matrix between equally-spaced grid cells on a torus, as some 
  # function `f` of euclidean distance. `x` and `y` are vectors containing 
  # equally-spaced vectors for the x and y coordinates of the grid cells on a
  # plane (i.e. the real coordinates). These should have been extended in order
  # to approximate some centre portion as a 2D plane.
  
  # number and dimension of grid cells
  m <- length(x)
  n <- length(y)
  x_size <- x[2] - x[1]
  y_size <- y[2] - y[1]
  
  # create indices for x and y on the first row of the dispersal matrix
  xidx <- rep(1:m, n)
  yidx <- rep(1:n, each = m)
  
  # project onto the torus and get distances from the first cell, in each
  # dimension
  xdist <- abs(xidx - xidx[1])
  ydist <- abs(yidx - yidx[1])
  xdist <- pmin(xdist, m - xdist) * x_size
  ydist <- pmin(ydist, n - ydist) * y_size
  
  # flatten distances into Euclidean space, apply the dispersal function and
  # return
  d <- sqrt(xdist ^ 2 + ydist ^ 2)
  f(d)
}


setupFFT <- function (x, y, f, factor = 2) {
  # set up the objects needed for the FFT dispersal matrix representation
  
  # extend the vectors (to project our plane on <= 1/4 of a torus)
  xe <- extend(x, factor)
  ye <- extend(y, factor)
  
  # get indices to true vectors
  xidx <- seq_range(attr(xe, 'idx'))
  yidx <- seq_range(attr(ye, 'idx'))
  
  # get fft basis for dispersal on a torus
  bcb_vec <- bcb(ye, xe, f)
  
  # create an empty population on all grid cells of the torus
  pop_torus <- matrix(0, length(ye), length(xe))
  
  # return as a named list for use in each iteration
  list(bcb_vec = bcb_vec,
       pop_torus = pop_torus,
       xidx = xidx,
       yidx = yidx)
} 


dispersalFFT <- function (popmat, fs) {
  # multiply the population matrix `popmat` giving the population of this stage 
  # in each cell through the dispersal matrix over the landscape, efficiently, 
  # and without ever having to construct the full matrix, by representing the 
  # landscape efficiently as a section of a torus and using linear algebra
  # identities of the fast Fourier transform
  
  # insert population matrix into the torus population
  fs$pop_torus[fs$yidx, fs$xidx] <- popmat
  
  # project population dispersal on the torus by fft
  
  # get spectral representations of the matrix & vector & compute the spectral
  # representation of their product
  pop_fft <- stats::fft(t(fs$pop_torus))
  bcb_fft <- stats::fft(fs$bcb_vec)
  pop_new_torus_fft <- ifft(pop_fft * bcb_fft)
  
  # convert back to real domain, apply correction and transpose
  pop_torus_new <- t(Re(pop_new_torus_fft / length(fs$pop_torus)))
  
  # extract the section of the torus representing our 2D plane and return
  pop_new <- pop_torus_new[fs$yidx, fs$xidx]
  pop_new
}


seq_range <- function (range, by = 1) seq(range[1], range[2], by = by)


ifft <- function (z) stats::fft(z, inverse = TRUE)


dispersal <- function(params, pop, hsm, cc, method){
  #stopifnot(is.dispersal(params))
  if(!any(method==c('ca','fft')))stop('method must be either "ca" (cellular automata) or\n "fft" (fast fourier transformation).')
  #stopifnot(is.habitat(habitat))  
  dispersal_results <- switch(method,
                              ca = dispersal_core_ca(params,pop,hsm,cc),
                              fft = dispersal_core_fft(params,pop))  
  return(dispersal_results)
}


dispersal_core_ca <- function(params, pop, hsm, cc){
  
  #generate default parameters for dispersal parameters if they are missing from 'params'. 
  if(!exists('barrier_type',params)) params$barrier_type <- 0
  if(!exists('dispersal_steps',params)) params$dispersal_steps <- 1
  if(!exists('use_barriers',params)) params$use_barriers <- FALSE
  
  #identify populations and workout which populations can disperse.
  which_stages_disperse <- which(params$dispersal_proportion>0)
  n_dispersing_stages <- length(which_stages_disperse)

  #if barriers is NULL create a barriers matrix all == 0.
  if(!exists('barriers_map',params)){
    bm <- raster::calc(hsm,function(x){x[!is.na(x)] <- 0; return(x)})
    params$barriers_map <- bm
  }
  
  # if(inherits(params$barriers_map,c("RasterStack","RasterBrick"))){
  #   bm <- params$barriers_map[[time_step]]
  #   params$barriers_map <- bm
  # }

  # could do this in parallel if wanted. 
  for (i in which_stages_disperse){
    pop[[i]][] <- rcpp_dispersal(raster::as.matrix(pop[[i]]),
                                      raster::as.matrix(cc),
                                      raster::as.matrix(hsm),
                                      raster::as.matrix(params$barriers_map),
                                      as.integer(params$barrier_type),
                                      params$use_barrier,
                                      as.integer(params$dispersal_steps),
                                      as.integer(params$dispersal_distance[i]),
                                      as.numeric(unlist(params$dispersal_kernel[i])),
                                      as.numeric(params$dispersal_proportion[i])
                                      )$dispersed_population
    
}

  return(pop)
}


dispersal_core_fft <- function(params, pop){
  
  #identify populations and workout which populations can disperse.
  which_stages_disperse <- which(params$dispersal_proportion>0)
  n_dispersing_stages <- length(which_stages_disperse)
  
  ## get the relevant 
  pops <- pop
  disperse_pops <- pops[which_stages_disperse]
  n <- dim(pops[[1]])[1:2]
  y <- seq_len(n[1])
  x <- seq_len(n[2])
  
  ## set up disperal function
  f <- function (d, cutoff = min(n)) {
    disp <- ifelse (d > cutoff, 0, exp(-d))
    disp / sum(disp)
  }
  
  # f <- function (d) exp(-d)
  # setup for the fft approach (run this once, before the simulation)
  fs <- setupFFT(x = x, y = y, f = f)
  
  # apply dispersal to the population (need to run this separately for each stage)
  # fft_dispersal <- list()
  # # could do this in parallel if wanted. 
  # for (i in seq_len(n_dispersing_stages)){
  #   fft_dispersal[[i]] <- dispersalFFT(popmat = raster::as.matrix(pops[[i]]), fs = fs)
  # }
  
  for (i in which_stages_disperse){
    pops[[i]][] <- dispersalFFT(popmat = raster::as.matrix(pops[[i]]), fs = fs)
  }
  
  # fft_dispersal <- lapply(fft_dispersal,function(x){pops[[1]][]<-x;return(pops[[1]])})
  # pops[[which_stages_disperse]] <- fft_dispersal
  # pops <- lapply(pops, `attr<-`, "habitat", "populations")
  return(pops)
}

####################################
### pre-defined module functions ###
####################################

#' @rdname population_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the fast_population_dynamics object to modify the  
#' # population using life-stage transitions and dispersal:
#' 
#' test_state2 <- fast_population_dynamics(test_state, 1)
#' 
#' par(mfrow=c(1,2))
#' plot(test_state$population$population_raster[[2]])
#' plot(test_state2$population$population_raster[[2]])

fast_population_dynamics <- function (state, timestep) {
  
  population_raster <- state$population$population_raster
  dispersal_parameters <- state$demography$dispersal_parameters
  transition_matrix <- state$demography$global_transition_matrix
  
  # get population as a matrix
  idx <- which(!is.na(raster::getValues(population_raster[[1]])))
  population <- raster::extract(population_raster, idx)
  
  # do population change
  population <- population %*% transition_matrix
  
  # do dispersal
  locations <- raster::xyFromCell(population_raster, idx)
  resolution <- mean(raster::res(population_raster))
  dispersal_decay <- dispersal_parameters * resolution
  
  dispersal <- dispersal_matrix(locations, dispersal_decay)
  population <- dispersal %*% population
  
  # put back in the raster
  population_raster[idx] <- population
  
  state$population$population_raster <- population_raster
  state
}

#' @rdname population_dynamics
#' 
#' @export
#' 
#' @examples
#'
#' # Use the ca_population_dynamics object to modify the  
#' # population using life-stage transitions, density-dependence,
#' # and cellular-automata based dispersal:
#' 
#' test_state2 <- ca_dispersal_population_dynamics(test_state_dp, 1)
#' 
#' par(mfrow=c(1,2))
#' plot(test_state_dp$population$population_raster[[2]])
#' plot(test_state_dp2$population$population_raster[[2]])

ca_dispersal_population_dynamics <- function (state, timestep) {
  
  population_raster <- state$population$population_raster
  dispersal_parameters <- state$demography$dispersal_parameters
  transition_matrix <- state$demography$global_transition_matrix
  transition_matrix_sd <- state$demography$transition_matrix_sd
  habitat_suitability <- state$habitat$habitat_suitability
  carrying_capacity <- state$habitat$carrying_capacity
  
  # get population as a matrix
  idx <- which(!is.na(raster::getValues(population_raster[[1]])))
  population <- raster::extract(population_raster, idx)
  
  # do population change
  population <- population %*% transition_matrix
  
  # check density dependence
  population <- cap_population(population, carrying_capacity)
  
  # put back in the raster
  population_raster[idx] <- population
   
  # do dispersal
  state$population$population_raster <- dispersal(params = dispersal_parameters,
                                                  pop = population_raster,
                                                  hsm = habitat_suitability,
                                                  cc = carrying_capacity,
                                                  method = "ca"
                                                  )
  state
} 

#' @rdname population_dynamics
#' 
#' @export
#'
#' @examples
#' 
#' test_state_dp2 <- fft_dispersal_population_dynamics(test_state_dp, 1)
#' 
#' par(mfrow=c(1,2))
#' plot(test_state_dp$population$population_raster[[2]])
#' plot(test_state_dp2$population$population_raster[[2]])

fft_dispersal_population_dynamics <- function (state, timestep) {
  
  population_raster <- state$population$population_raster
  dispersal_parameters <- state$demography$dispersal_parameters
  transition_matrix <- state$demography$global_transition_matrix
  habitat_suitability <- state$habitat$habitat_suitability
  carrying_capacity <- state$habitat$carrying_capacity
  
  # get population as a matrix
  idx <- which(!is.na(raster::getValues(population_raster[[1]])))
  population <- raster::extract(population_raster, idx)
  
  # do population change
  population <- population %*% transition_matrix
  
  # put back in the raster
  population_raster[idx] <- population
  
  # do dispersal
  state$population$population_raster <- dispersal(params = dispersal_parameters,
                                                  pop = population_raster,
                                                  method = "fft"
  )
  state
}
