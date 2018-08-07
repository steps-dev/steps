#' @useDynLib steps
#' @importFrom Rcpp sourceCpp
NULL
#' Functions to modify the population in a state object.
#' 
#' Pre-defined functions to operate on a population
#' during a simulation.
#'
#' @rdname population_dynamics_functions
#'
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
#' @examples
#' 
#' library(steps)
#' 
#' @rdname population_dynamics_functions
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
        
        local_t <- array(apply(local_mat, 3, function(x) rbind(0, x[-1,], 1 - apply(x[-1, ], 2, sum))),
                         dim = c(dim(local_mat)[1]+1, dim(local_mat)[2], dim(local_mat)[3]))
        local_f <- local_mat
        local_f[-1, , ] <- 0
        
        pop_tmp <- cbind(population, rep(0, nrow(population)))
        
        survival_stochastic <- sapply(seq_len(ncol(population)),
                                      function(y) sapply(seq_len(nrow(population)),
                                                         function(x) rmultinom(n = 1,
                                                                            size = pop_tmp[x, y],
                                                                            prob = local_t[, y, x])),
                                      simplify = 'array')
                                      

        new_offspring_deterministic <- sapply(seq_len(nrow(population)), function(x) local_f[ , , x] %*% matrix(population[x, ]))
        new_offspring_stochastic <- matrix(rpois(n = length(c(new_offspring_deterministic)),
                                                 lambda = c(new_offspring_deterministic)),
                                           nrow = nrow(new_offspring_deterministic))
        new_offspring <- apply(new_offspring_stochastic, 2, sum)

        population <- t(apply(survival_stochastic[seq_len(ncol(population)), , ], c(1, 2), sum))
        population[ , 1] <- population[ , 1] + new_offspring
        
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
                             function(x) rmultinom(n = nrow(population),
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

#' @rdname population_dynamics_functions
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


#' @rdname population_dynamics_functions
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

kernel_dispersal <- function(
  dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.1),
  dispersal_proportion = list(0, 0.35, 0.35 * 0.714, 0),
  arrival_probability = "both",
  stages = NULL,
  fft = FALSE
) {
  
  if (fft & !is.null(arrival_probability)) {
    stop(
      "An arrival probability surface can't be used with FFT approximation.\n",
      "Set \"arrival_probability = NULL\""
    )
  }

  pop_dynamics <- function(state, timestep) {

    # Which stages can disperse
    which_stages_disperse <- which(dispersal_proportion > 0)
    
    # Which stages contribute to density dependence.
    which_stages_density <- if (is.null(stages)) {
        seq(raster::nlayers(state$population$population_raster))
      } else {
        stages
      }
    
    if (fft) {
      # Apply dispersal to the population
      # (need to run this separately for each stage)
      for (stage in which_stages_disperse) {
        state$population$population_raster[[stage]][] <- dispersalFFT(
          popmat = raster::as.matrix(
            state$population$population_raster[[stage]]
          ),
          fs = setupFFT(
            x = raster::ncol(state$population$population_raster),
            y = raster::nrow(state$population$population_raster),
            f = function(d) {
                  disp <- dispersal_kernel(d)
                  disp / sum(disp)
                }
          )
        )
      }

    } else {

      # Extract locations as x and y coordinates from landscape (ncells x 2)
      xy <- raster::xyFromCell(
        state$population$population_raster,
        seq(raster::ncell(state$population$population_raster))
      )

      # Rescale locations
      xy <- sweep(xy, 2, raster::res(state$population$population_raster), "/")

      # Extract arrival probabilities
      arrival_probability <- match.arg(
        arrival_probability,
        c("both", "habitat_suitability", "carrying_capacity")
      )

      delayedAssign(
        "habitat_suitability_values",
        raster::getValues(state$habitat$habitat_suitability)
      )

      delayedAssign(
        "carrying_capacity_proportion",
        raster::getValues(
          raster::calc(
            raster::stack(state$population$population_raster)[[
              which_stages_density
            ]],
            sum
          ) /
          state$habitat$carrying_capacity
        )
      )

      arrival_prob_values <- switch(
        arrival_probability,
        both = habitat_suitability_values * carrying_capacity_proportion,
        habitat_suitability = habitat_suitability_values,
        carrying_capacity = carrying_capacity_proportion
      )
      
      # Only non-zero arrival prob cells can receive individuals
      can_arriv <- which(arrival_prob_values > 0 & !is.na(arrival_prob_values))

      for (stage in which_stages_disperse) {

        # Extract the population values
        population_values <- raster::getValues(
          state$population$population_raster[[stage]]
        )

        # Only non-zero population cells can contribute
        has_pop <- which(population_values > 0 & !is.na(population_values))

        contribute <- function(i) {
          # Euclidean distance between ith cell and all the cells that it can
          # contribute to
          contribution <- sqrt(
            (xy[i, "x"] - xy[can_arriv, "x"])^2 +
            (xy[i, "y"] - xy[can_arriv, "y"])^2
          )
          contribution <- dispersal_kernel(contribution)
          contribution <- contribution * arrival_prob_values[can_arriv]
          
          # Standardise contributions
          contribution <- contribution / sum(contribution)
          contribution * population_values[i]
        }

        state$population$population_raster[[stage]][can_arriv] <- rowSums(
          vapply(has_pop, contribute, as.numeric(can_arriv))
        )
      }
    }

    state

}

  as.population_kernel_dispersal(pop_dynamics)

}


#' @rdname population_dynamics_functions
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

  
#' @rdname population_dynamics_functions
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


#' @rdname population_dynamics_functions
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