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
#' @param pop_change a function to define how population growth occurs at each timestep
#' @param pop_disp a function to define how the population disperses at each timestep (default is exponential kernel)
#' @param pop_mod a function to define any deterministic changes to the population - such as translocation - at each timestep
#' @param pop_dens_dep a function to control density dependence effects on the population at each timestep
#' @param demo_stoch should demographic stochasticity be used in population change? (default is FALSE)
#' @param dispersal_kernel a single or list of user-defined distance dispersal kernel functions
#' @param dispersal_proportion proportions of individuals (0 to 1) that can disperse in each life stage
#' @param arrival_probability a raster layer that controls where individuals can disperse to (e.g. habitat suitability)
#' @param fft should a fast-fourier approximation be used?
#' @param dispersal_distance the distances (in cell units) that each life stage can disperse
#' @param barrier_type if barrier map is used, does it stop (0 - default) or kill (1) individuals
#' @param dispersal_steps number of dispersal steps to take before stopping
#' @param use_barriers should dispersal barriers be used? If so, a barriers map must be provided
#' @param barriers_map a raster layer that contains cell values of 0 (no barrier) and 1 (barrier)
#' @param carrying_capacity a raster layer that specifies the carrying capacity in each cell
#' @param source_layer a spatial layer with the locations and number of individuals to translocate from - note, this layer will only have zero values if individuals are being introduced from outside the study area
#' @param sink_layer a spatial layer with the locations and number of individuals to translocate to
#' @param effect_timesteps which timesteps in a single simulation do the translocations take place
#' @param stages which life-stages contribute to density dependence or are affected by the translocations - default is all
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
#' pop_source <- r
#' pop_source[!is.na(r)] <- 0
#' pop_source[!is.na(r) & r >= 1200] <- 100
#' 
#' pop_sink <- r
#' pop_sink[!is.na(r)] <- 0
#' pop_sink[sample(which(!is.na(getValues(r)) & getValues(r) < 200),
#'                 length(pop_source[!is.na(r) & r >= 1200]))] <- 100
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
#' test_demography <- build_demography(transition_matrix = mat)
#' 
#' test_population <- build_population(pop)
#' 
#' test_state <- build_state(test_habitat, test_demography, test_population)
#'                              
#' example_function <- function (state, timestep) {
#'   state
#' }
#' 
#' example_function <- as.population_dynamics(example_function)

as.population_dynamics <- function (population_dynamics_function) {
  as_class(population_dynamics_function, "population_dynamics", "function")
}


#' @rdname population_dynamics
#'
#' @export
#' 
#' @examples
#'
#' # Test if object is of the type 'population dynamics'
#'   
#' is.population_dynamics(example_function)

is.population_dynamics <- function (x) {
  inherits(x, 'population_dynamics')
}

#' @rdname population_dynamics
#'
#' @export
#'
#' @examples
#' 
#' print(example_function)

print.population_dynamics <- function (x, ...) {
  cat("This is a population_dynamics object")
}


#' @rdname population_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the population_dynamics object to modify the  
#' # population with population change, dispersal, density dependence,
#' # and population modification functions:
#'
#' example_function <- population_dynamics()
#' test_state2 <- example_function(test_state, 1)
#' 
#' par(mfrow=c(1,2))
#' plot(test_state$population$population_raster[[2]])
#' plot(test_state2$population$population_raster[[2]])

population_dynamics <- function (pop_change = simple_growth(),
                                 pop_disp = NULL,
                                 pop_mod = NULL,
                                 pop_dens_dep = NULL) {
  
  pop_dynamics <- function (state, timestep) {
    
    state <- pop_change(state, timestep)
    
    if (!is.null(pop_disp))
      state <- pop_disp(state, timestep)
    
    if (!is.null(pop_mod))
      state <- pop_mod(state, timestep)
    
    if (!is.null(pop_dens_dep))
      state <- pop_dens_dep(state, timestep)

    state
  }
  
  as.population_dynamics(pop_dynamics)
  
}

##########################
### internal functions ###
##########################

as.population_simple_growth <- function (population_simple_growth) {
  as_class(population_simple_growth, "population_dynamics", "function")
}

as.population_demo_stoch <- function (population_demo_stoch) {
  as_class(population_demo_stoch, "population_dynamics", "function")
}

as.population_kernel_dispersal <- function (population_kernel_dispersal) {
  as_class(population_kernel_dispersal, "population_dynamics", "function")
}

as.population_ca_dispersal <- function (population_ca_dispersal) {
  as_class(population_ca_dispersal, "population_dynamics", "function")
}

as.population_fft_dispersal <- function (population_fft_dispersal) {
  as_class(population_fft_dispersal, "population_dynamics", "function")
}

as.population_translocation <- function (population_translocation) {
  as_class(population_translocation, "population_dynamics", "function")
}

as.population_density_dependence <- function (population_density_dependence) {
  as_class(population_density_dependence, "population_dynamics", "function")
}


###### DISPERSAL FUNCTIONS ######

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
  n2 <- 2 ^ ceiling(log2(factor * n))
  
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

  pop_new_torus_fft <- stats::fft(pop_fft * bcb_fft, inverse = TRUE)
  
  # convert back to real domain, apply correction and transpose
  pop_torus_new <- t(Re(pop_new_torus_fft) / length(fs$pop_torus))
  pop_torus_new <- pmax(pop_torus_new, 0)
  
  # extract the section of the torus representing our 2D plane and return
  pop_new <- pop_torus_new[fs$yidx, fs$xidx]

  # check for missing values and replace with zeros
  if (any(is.na(pop_new))) {
    pop_new[is.na(pop_new)] <- 0
  }
  
  # make sure none are lost or gained (unless all are zeros)
  if (any(pop_new[] > 0)) {
    pop_new[] <- stats::rmultinom(1, size = sum(popmat), prob = pop_new[])    
  }
  
  pop_new
}


seq_range <- function (range, by = 1) seq(range[1], range[2], by = by)

# dispersal <- function(params, pop, hsm, cc, method){
#   #stopifnot(is.dispersal(params))
#   if(!any(method==c('ca','fft')))stop('method must be either "ca" (cellular automata) or\n "fft" (fast fourier transformation).')
#   #stopifnot(is.habitat(habitat))  
#   dispersal_results <- switch(method,
#                               ca = dispersal_core_ca(params,pop,hsm,cc),
#                               fft = dispersal_core_fft(params,pop))  
#   return(dispersal_results)
# }


# dispersal_core_ca <- function(dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
#                               dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
#                               dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
#                               barrier_type = 0,
#                               dispersal_steps = 1,
#                               use_barriers = FALSE,
#                               barriers_map = NULL,
#                               arrival_probability = "habitat_suitability",
#                               carrying_capacity = "carrying_capacity",
#                               pop){
#   
#   #generate default parameters for dispersal parameters if they are missing from 'params'. 
#   #if(!exists('barrier_type',params)) barrier_type <- 0
#   #if(!exists('dispersal_steps',params)) params$dispersal_steps <- 1
#   #if(!exists('use_barriers',params)) params$use_barriers <- FALSE
#   
#   #identify populations and workout which populations can disperse.
#   which_stages_disperse <- which(dispersal_proportion>0)
#   n_dispersing_stages <- length(which_stages_disperse)
# 
#   #if barriers is NULL create a barriers matrix all == 0.
#   if(is.null(barriers_map)){
#     bm <- raster::calc(state$habitat[[arrival_probability]],
#                        function(x){x[!is.na(x)] <- 0; return(x)})
#     barriers_map <- bm
#   }
#   
#   # if(inherits(params$barriers_map,c("RasterStack","RasterBrick"))){
#   #   bm <- params$barriers_map[[time_step]]
#   #   params$barriers_map <- bm
#   # }
# 
#   # could do this in parallel if wanted. 
#   for (i in which_stages_disperse){
#     pop[[i]][] <- rcpp_dispersal(raster::as.matrix(pop[[i]]),
#                                       raster::as.matrix(state$habitat[[carrying_capacity]]),
#                                       raster::as.matrix(state$habitat[[arrival_probability]]),
#                                       raster::as.matrix(barriers_map),
#                                       as.integer(barrier_type),
#                                       use_barrier,
#                                       as.integer(dispersal_steps),
#                                       as.integer(dispersal_distance[i]),
#                                       as.numeric(unlist(dispersal_kernel[i])),
#                                       as.numeric(dispersal_proportion[i])
#                                       )$dispersed_population
#     
# }
# 
#   return(pop)
# }


# dispersal_core_fft <- function(dispersal_proportion,
#                                pop){
#   
#   #identify populations and workout which populations can disperse.
#   which_stages_disperse <- which(dispersal_proportion>0)
#   n_dispersing_stages <- length(which_stages_disperse)
# 
#   ## get the relevant 
#   pops <- pop
#   disperse_pops <- pops[which_stages_disperse]
#   n <- dim(pops[[1]])[1:2]
#   y <- seq_len(n[1])
#   x <- seq_len(n[2])
#   
#   ## set up disperal function
#   f <- function (d, max = Inf, mean = 1) {
#     lambda <- 1 / mean
#     disp <- ifelse(d > max,
#                    0,
#                    exp(-lambda * d))
#     disp / sum(disp)
#   }
# 
#   # setup for the fft approach (run this once, before the simulation)
#   fs <- setupFFT(x = x, y = y, f = f)
#   
#   # apply dispersal to the population (need to run this separately for each stage)
# 
#   for (i in which_stages_disperse){
#     pops[[i]][] <- dispersalFFT(popmat = raster::as.matrix(pops[[i]]), fs = fs)
#   }
# 
#   return(pops)
# }

####################################
### pre-defined module functions ###
####################################

#' @rdname population_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the simple growth function to modify the  
#' # population using life-stage transitions:
#'
#' test_lin_growth <- simple_growth()

simple_growth <- function (demo_stoch = FALSE) {
  
  pop_dynamics <- function (state, timestep) {
    
    # import components from state object
    population_raster <- state$population$population_raster
    transition_matrix <- state$demography$global_transition_matrix
    #nstages <- ncol(state$demography$global_transition_matrix) not used...

    # get population as a matrix
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population <- raster::extract(population_raster, idx)
    
    # do population change
    if (!is.null(state$demography$local_transition_matrix)) {
      
      local_mat <-  state$demography$local_transition_matrix
      
      if(demo_stoch){
        
        
        
      } else {
        
        population <- t(sapply(seq_len(nrow(population)), function(x) local_mat[ , , x] %*% matrix(population[x, ])))
        
      }

    } else {
      #population <- apply(population, 1, function(x) transition_matrix %*% as.matrix(x))
      
      if(demo_stoch){
        
        t <- rbind(0, transition_matrix[-1,], 1 - apply(transition_matrix[-1, ], 2, sum))
        f <- transition_matrix
        f[-1, ] <- 0
        
        # population <- sapply(seq_len(ncol(population)),
        #                             function(x) stats::rbinom(nrow(population),
        #                                                       population[, x],
        #                                                       t[which(t[, x] > 0), x])
        #                             )

        pop_tmp <- cbind(population, rep(0, nrow(population)))
        
        survival_stochastic <- sapply(seq_len(ncol(population)),
                             function(x) stats::rmultinom(n = nrow(population),
                                                       size = pop_tmp[, x],
                                                       prob = t[, x]),
                             simplify = 'array'
        )
        
        new_offspring_deterministic <- f %*% t(population)
        new_offspring_stochastic <- matrix(rpois(n = length(c(new_offspring_deterministic)),
                                          lambda = c(new_offspring_deterministic)),
                                          nrow = nrow(new_offspring_deterministic))
        new_offspring <- apply(new_offspring_stochastic, 2, sum)
        
        population <- t(apply(survival_stochastic[seq_len(ncol(population)), , ], c(1, 2), sum))
        population[ , 1] <- population[ , 1] + new_offspring
        
        

        # #fecundity
        # newborns <- stats::rpois(n,
        #                          transition_matrix[1, j] * population[ , j]) 
        # #survival
        # survivors <- stats::rbinom(n,
        #                            ceiling(as.vector(population[ , j])),
        #                            transition_matrix[utils::tail(which(transition_matrix[ , j] > 0),1) , j])
        
      } else {
        
        population <- t(transition_matrix %*% t(population))

      }

    }

    # put back in the raster
    population_raster[idx] <- population
    
    state$population$population_raster <- population_raster
    
    state
  }
  
  as.population_simple_growth(pop_dynamics)
  
}

#' @rdname population_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the demographic stochasticity function to modify the  
#' # population using random variation:
#'
#' test_dem_stoch <- demographic_stochasticity()

demographic_stochasticity <- function () {
  
  pop_dynamics <- function (state, timestep) {
    
    population_raster <- state$population$population_raster
    transition_matrix <- state$demography$global_transition_matrix
    
    # get population as a matrix
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population <- raster::extract(population_raster, idx)

    if (!is.null(state$demography$local_transition_matrix)) {
      
      local_mat <-  state$demography$local_transition_matrix
      
      for (i in seq_len(nrow(population))) {
        for (j in seq_len(ncol(population))) {
          
          #fecundity
          newborns <- stats::rpois(1, local_mat[1, j, i] * population[i , j]) 
          
          #survival
          survivors <- stats::rbinom(1,
                                     ceiling(as.vector(population[i , j])),
                                     local_mat[utils::tail(which(local_mat[ ,j , i] > 0),1) , j, i])
          
          population[i, j] <- newborns + survivors
          
        }
      }

    } else {
    
      n <- nrow(population)
      
      for (j in seq_len(ncol(population))) {
        
        #fecundity
        newborns <- stats::rpois(n,
                                 transition_matrix[1, j] * population[ , j]) 
        #survival
        survivors <- stats::rbinom(n,
                                   ceiling(as.vector(population[ , j])),
                                   transition_matrix[utils::tail(which(transition_matrix[ , j] > 0),1) , j])
        
        #population[, j] <- newborns + survivors
        population[ , j] <- survivors
        population[ , 1] <- newborns
        
      }
    }
    
    population_raster[idx] <- population
    
    state$population$population_raster <- population_raster
    
    state
  }
  
  as.population_demo_stoch(pop_dynamics)
  
}


# #' @rdname population_dynamics
# #' 
# #' @export
# #' 
# #' @examples
# #' 
# #' # Use the simple dispersal function to modify the  
# #' # population using a diffusion kernel:
# #'
# #' test_sim_dispersal <- simple_dispersal(dispersal_parameters = dp_params)

# simple_dispersal <- function (distance_decay = 0.5, dispersal_parameters) {
# 
#   pop_dynamics <- function (state, timestep) {
# 
#     # import components from state object
#     population_raster <- state$population$population_raster
# 
#     # identify populations and workout which populations can disperse.
#     which_stages_disperse <- which(dispersal_parameters$dispersal_proportion>0)
#     n_dispersing_stages <- length(which_stages_disperse)
#     
#     # get population as a matrix - one column per life stage, one row per grid cell
#     idx <- which(!is.na(raster::getValues(population_raster[[1]])))
#     population <- raster::extract(population_raster, idx)[, which_stages_disperse]
# 
#     # extract locations as x and y coordinates from landscape
#     locations <- raster::xyFromCell(population_raster, idx)
#     
#     # get resolution and rescale distance decay accordingly
#     resolution <- mean(raster::res(population_raster))
#     dispersal_decay <- distance_decay * resolution
# 
#     # calculate distance matrix
#     D <- as.matrix(stats::dist(locations))
#     
#     # redistribute populations based on probability density function
#     dispersal_matrix <- exp(-D / dispersal_decay)
#     rm(D)
#     
#     # get population sums in each cell and rescale to maintain total populations
#     sums <- colSums(dispersal_matrix)
#     dispersal <- sweep(dispersal_matrix, 2, sums, "/")
#     rm(dispersal_matrix)
#     
#     # multiply populations by new dispersal locations
#     population <- dispersal %*% population
# 
#     # put back in the raster
#     for (i in seq_len(n_dispersing_stages)){
#       population_raster[[as.integer(which_stages_disperse)[i]]][idx] <- population[ , i]
#     }
#     #population_raster[[which_stages_disperse]][idx] <- population
#     state$population$population_raster <- population_raster
#     
#     state
#   }
# 
#   as.population_simple_dispersal(pop_dynamics)
# 
# }


#' @rdname population_dynamics
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the kernel-based dispersal function to modify the  
#' # population using a user-defined diffusion distribution
#' # and an optional arrival probability layer (e.g. habitat suitability):
#'
#' test_kern_dispersal <- kernel_dispersal()

kernel_dispersal <- function (dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.1),
                              dispersal_proportion = list(0, 0.35, 0.35*0.714, 0),
                              arrival_probability = "habitat_suitability",
                              fft = FALSE
                              ) {
  
  if (fft & !is.null(arrival_probability)) {
    stop("An arrival probability surface can't be used with the FFT approximation.\n",
         "Set arrival_probability = NULL in the kernel_dispersal function parameters")
  }
  
  pop_dynamics <- function (state, timestep) {
    
    # import components from state object
    population_raster <- state$population$population_raster
    
    if (!is.null(arrival_probability)) {
      arrival_prob <- state$habitat[[arrival_probability]]      
    }

    # identify populations and workout which populations can disperse
    which_stages_disperse <- which(dispersal_proportion > 0)
    n_dispersing_stages <- length(which_stages_disperse)
    
    # get population as a matrix - one column per life stage, one row per grid cell
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population <- raster::extract(population_raster, idx)
    
    # get cell arrival probabilities as a vector
    if (!is.null(arrival_probability)) {
      suitabilities <- raster::extract(arrival_prob, idx)
    }

    if(fft) {

      disperse_pops <- population_raster[which_stages_disperse]
      n <- dim(population_raster[[1]])[1:2]
      y <- seq_len(n[1])
      x <- seq_len(n[2])
      
      ## set up disperal function
      f <- function (d) {
        disp <- dispersal_kernel(d)
        disp / sum(disp)
      }
      
      # setup for the fft approach (run this once, before the simulation)
      fs <- setupFFT(x = x, y = y, f = f)
      
      # apply dispersal to the population (need to run this separately for each stage)
      for (i in which_stages_disperse){
        population_raster[[i]][] <- dispersalFFT(popmat = raster::as.matrix(population_raster[[i]]), fs = fs)
      }
      
      state$population$population_raster <- population_raster
      
    } else {
      
      # extract locations as x and y coordinates from landscape (ncells x 2)
      locations <- raster::xyFromCell(population_raster, idx)
      
      # determine the cell resolution for calculating cell distances traveled
      resolution <- mean(raster::res(population_raster))
      
      # construct a matrix all distances between all cells (ncells x ncells)
      # rescale based on cell resolution
      D <- as.matrix(stats::dist(locations))/resolution
      
      # apply dispersal function to determine proportions of individuals that remain/move into cells
      # adjust the distance by the raster resolution to convert to cell distances
      dispersal_matrix <- dispersal_kernel(D)
      
      dispersal_matrix <- sweep(dispersal_matrix, 1, suitabilities, "*")
      
      # rescale column values to sum to one
      sums <- colSums(dispersal_matrix)
      dispersal <- sweep(dispersal_matrix, 2, sums, "/")
      
      # determine dispersing population values
      population <- dispersal %*% population
      
      # put back in the raster
      for (i in seq_len(n_dispersing_stages)){
        population_raster[[as.integer(which_stages_disperse)[i]]][idx] <- population[ , i]
      }
      state$population$population_raster <- population_raster
      
    }
 
    state
  }
  
  as.population_kernel_dispersal(pop_dynamics)
  
}


#' @rdname population_dynamics
#' 
#' @export
#' 
#' @examples
#'
#' # Use the cellular automata dispersal function to modify  
#' # the population using rule-based cell movements:
#' 
#' test_ca_dispersal <- cellular_automata_dispersal()

cellular_automata_dispersal <- function (dispersal_distance=list(0, 10, 10, 0),
                                         dispersal_kernel=list(0, exp(-c(0:9)^1/3.36), exp(-c(0:9)^1/3.36), 0),
                                         dispersal_proportion=list(0, 0.35, 0.35*0.714, 0),
                                         barrier_type = 0,
                                         dispersal_steps = 1,
                                         use_barriers = FALSE,
                                         barriers_map = NULL,
                                         arrival_probability = "habitat_suitability",
                                         carrying_capacity = "carrying_capacity") {

  pop_dynamics <- function (state, timestep) {

    population_raster <- state$population$population_raster
    arrival_prob <- state$habitat[[arrival_probability]]
    carrying_capacity <- state$habitat[[carrying_capacity]]

    # get population as a matrix
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population <- raster::extract(population_raster, idx)

    # identify dispersing stages
    which_stages_disperse <- which(dispersal_proportion>0)
    n_dispersing_stages <- length(which_stages_disperse)
    
    #if barriers is NULL create a barriers matrix all == 0.
    if(is.null(barriers_map)){
      barriers_map <- raster::calc(arrival_prob,
                         function(x){x[!is.na(x)] <- 0; return(x)})
    }
    
    # if(inherits(params$barriers_map,c("RasterStack","RasterBrick"))){
    #   bm <- params$barriers_map[[time_step]]
    #   params$barriers_map <- bm
    # }
    
    # could do this in parallel if wanted. 
    for (i in which_stages_disperse){
      population_raster[[i]][] <- rcpp_dispersal(raster::as.matrix(population_raster[[i]]),
                                   raster::as.matrix(carrying_capacity),
                                   raster::as.matrix(arrival_prob),
                                   raster::as.matrix(barriers_map),
                                   as.integer(barrier_type),
                                   use_barriers,
                                   as.integer(dispersal_steps),
                                   as.integer(dispersal_distance[i]),
                                   as.numeric(unlist(dispersal_kernel[i])),
                                   as.numeric(dispersal_proportion[i])
      )$dispersed_population
      
    }
 
    state$population$population_raster <- population_raster
    
    state
  }

  as.population_ca_dispersal(pop_dynamics)

}

  
# #' @rdname population_dynamics
# #' 
# #' @export
# #'
# #' @examples
# #'
# #' # Use the fast fourier dispersal function to modify the  
# #' # population using rule-based cell movements:
# #' 
# #' test_fft_dispersal <- fast_fourier_dispersal(dispersal_parameters = dp_params)

# fast_fourier_dispersal <- function (dispersal_parameters) {
# 
#   pop_dynamics <- function (state, timestep) {
# 
#     population_raster <- state$population$population_raster
#     habitat_suitability <- state$habitat$habitat_suitability
#     carrying_capacity <- state$habitat$carrying_capacity
# 
#     # get population as a matrix
#     idx <- which(!is.na(raster::getValues(population_raster[[1]])))
#     population <- raster::extract(population_raster, idx)
# 
#     # do dispersal
#     state$population$population_raster <- dispersal_core_fft(params = dispersal_parameters,
#                                                     pop = population_raster)
#     state
#   }
# 
#   as.population_fft_dispersal(pop_dynamics)
# 
# }


#' @rdname population_dynamics
#'
#' @export
#' 
#' @examples
#' 
#' # Use the translocation_population_dynamics object to modify the  
#' # population using translocations:
#' 
#' test_ca_dispersal <- pop_translocation(source_layer = pop_source,
#'                                        sink_layer = pop_sink,
#'                                        stages = NULL,
#'                                        effect_timesteps = 1)

pop_translocation <- function (source_layer, sink_layer, stages = NULL, effect_timesteps = NULL) {
  
  pop_dynamics <- function (state, timestep) {
    
    if (timestep %in% effect_timesteps) {
      
      population_raster <- state$population$population_raster
      nstages <- ncol(state$demography$global_transition_matrix)
      
      # get population as a matrix
      idx <- which(!is.na(raster::getValues(population_raster[[1]])))
      population_matrix <- raster::extract(population_raster, idx)
      
      source <- raster::extract(source_layer, idx)
      sink <- raster::extract(sink_layer, idx)
      
      if (is.null(stages)) {
        
        for (stage in seq_len(ncol(population_matrix))) {

          if (any(population_matrix[ , stage] + (sink/nstages) - (source/nstages) < 0)) {
            stop("Your translocations are resulting in negative populations -\nsource locations/numbers are not viable.")
          }
          population_matrix[ , stage] <- population_matrix[ , stage] + (sink/nstages) - (source/nstages)

        }

      } else {
        
        for (stage in stages) {
          if (any(population_matrix[ , stage] + (sink/length(stages)) - (source/length(stages)) < 0)) {
            stop("Your translocations are resulting in negative populations -\nsource locations/numbers are not viable.")
          }
          population_matrix[ , stage] <- population_matrix[ , stage] + (sink/length(stages)) - (source/length(stages))

        }
        
      #population_matrix[population_matrix < 0] <- 0 
        
      }

      # put back in the raster
      population_raster[idx] <- population_matrix

      state$population$population_raster <- population_raster

    }
      
    state
    
  }
  
  as.population_translocation(pop_dynamics)
  
}


#' @rdname population_dynamics
#'
#' @export
#' 
#' @examples
#' 
#' # Use the translocation_population_dynamics object to modify the  
#' # population using translocations:
#' 
#' test_pop_dd <- pop_density_dependence()

pop_density_dependence <- function (stages = NULL) {
  
  pop_dynamics <- function (state, timestep) {
    
    population_raster <- state$population$population_raster
    carrying_capacity <- state$habitat$carrying_capacity

    # get population as a matrix
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population_matrix <- raster::extract(population_raster, idx)
     
    # if (is.null(carrying_capacity)) {
    #   stop ("carrying capacity must be specified",
    #         call. = FALSE)
    # }
    
    # get degree of overpopulation, and shrink accordingly
    if (!is.null(stages)) {
      overpopulation <- as.vector(carrying_capacity) / rowSums(population_matrix[ ,stages])
    } else {
      overpopulation <- as.vector(carrying_capacity) / rowSums(population_matrix)
    }
    
    overpopulation[is.nan(overpopulation)] <- 0
    overpopulation <- pmin(overpopulation, 1)
    population <- sweep(population_matrix, 1, overpopulation, "*")
    
    # put back in the raster
    population_raster[idx] <- population
    
    state$population$population_raster <- population_raster 
      
    state
  }

  as.population_density_dependence(pop_dynamics)
  
}
