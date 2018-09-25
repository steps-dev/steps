#' @useDynLib steps
#' @importFrom Rcpp sourceCpp
NULL
#' Functions to modify the population in a landscape object.
#'
#' Pre-defined functions to operate on a population during a simulation.
#'
#' @name population_dynamics_functions
#'
#' @param transition_matrix A symmetrical age-based (Leslie) or stage-based
#'   population structure matrix.
#' @param demographic_stochasticity should demographic stochasticity be used in
#'   population change?
#' @param global_stochasticity,local_stochasticity either scalar values or
#'   matrices (with the same dimension as \code{transition_matrix}) giving
#'   variability (in standard deviations) in the transition matrix either for
#'   all populations (\code{global_stochasticity}) or for each population
#'   separately (\code{local_stochasticity})
#' @param dispersal_kernel a single or list of user-defined distance dispersal
#'   kernel functions
#' @param dispersal_proportion proportions of individuals (0 to 1) that can
#'   disperse in each life stage
#' @param distance_function defines distance between source cells and all
#'   potential sink cells for dispersal
#' @param arrival_probability a raster layer that controls where individuals can
#'   disperse to (e.g. habitat suitability)
#' @param dispersal_distance the distances (in cell units) that each life stage
#'   can disperse
#' @param stages which life-stages disperse, are modified (e.g. translocated),
#'   or contribute to density dependence - default is all
#' @param barrier_type if barrier map is used, does it stop (0 - default) or
#'   kill (1) individuals
#' @param dispersal_steps number of dispersal steps to take before stopping
#' @param use_barriers should dispersal barriers be used? If so, a barriers map
#'   must be provided
#' @param barriers_map a raster layer that contains cell values of 0 (no
#'   barrier) and 1 (barrier)
#' @param carrying_capacity a raster layer that specifies the carrying capacity
#'   in each cell
#' @param source_layer a spatial layer with the locations and number of
#'   individuals to translocate from - note, this layer will only have zero
#'   values if individuals are being introduced from outside the study area
#' @param sink_layer a spatial layer with the locations and number of
#'   individuals to translocate to
#' @param effect_timesteps which timesteps in a single simulation do the
#'   translocations take place
#'
#' @examples
#'
#' library(steps)


#' @rdname population_dynamics_functions
#' 
#' @export
#' 
#' @examples
#' 
#' # Use a function to modify the  
#' # population using life-stage transitions:
#'
#' test_growth <- growth(egk_mat)
growth <- function (transition_matrix,
                    demographic_stochasticity = TRUE,
                    global_stochasticity = 0,
                    local_stochasticity = 0) {
  
  idx <- which(transition_matrix != 0)
  is_recruitment <- upper.tri(transition_matrix)[idx]
  upper <- ifelse(is_recruitment, Inf, 1)
  vals <- transition_matrix[idx]
  dim <- nrow(transition_matrix)
  
  if (is.matrix(global_stochasticity)) {
    stopifnot(identical(dim(transition_matrix), dim(global_stochasticity)))
    stopifnot(identical(which(global_stochasticity != 0), idx))
    global_stochasticity <- global_stochasticity[idx]
  }

  if (is.matrix(local_stochasticity)) {
    stopifnot(identical(dim(transition_matrix), dim(local_stochasticity)))
    stopifnot(identical(which(local_stochasticity != 0), idx))
    local_stochasticity <- local_stochasticity[idx]
  }
  
  pop_dynamics <- function (landscape, timestep) {
    
    # import components from landscape object
    population_raster <- landscape$population
    
    # get population as a matrix
    cell_idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population <- raster::extract(population_raster, cell_idx)
    
    n_cells <- length(cell_idx)
    
    # add global noise to the transition matrices and truncate
    global_noise <- stats::rnorm(length(idx), 0, global_stochasticity)
    local_noise <- stats::rnorm(length(idx) * n_cells, 0, local_stochasticity)
    total_noise <- global_noise + local_noise
    
    values <- vals + total_noise
    values <- pmax(values, 0)
    values <- pmin(values, rep(upper, n_cells))
    
    transition_array <- array(0, c(dim, dim, n_cells))
    
    # pad the index to get corresponding elements in each slice  
    addition <- length(transition_matrix) * (seq_len(n_cells) - 1)
    idx_full <- as.numeric(outer(idx, addition, FUN = "+"))
    transition_array[idx_full] <- values

    if (demographic_stochasticity) {
      
      local_t <- array(apply(transition_array, 3,
                             function(x) {
                               rbind(0,
                                     x[-1,],
                                     1 - apply(x[-1, ], 2, sum))
                             }),
                       dim = c(dim + 1,
                               dim,
                               n_cells))
      
      local_f <- transition_array
      local_f[-1, , ] <- 0
      
      pop_tmp <- cbind(population, rep(0, n_cells))
      
      survival_stochastic <- sapply(seq_len(ncol(population)),
                                    function(y) sapply(seq_len(nrow(population)),
                                                       function(x) stats::rmultinom(n = 1,
                                                                                    size = pop_tmp[x, y],
                                                                                    prob = local_t[, y, x])),
                                    simplify = 'array')
      
      
      new_offspring_deterministic <- sapply(seq_len(nrow(population)), function(x) local_f[ , , x] %*% matrix(population[x, ]))
      new_offspring_stochastic <- matrix(stats::rpois(n = length(c(new_offspring_deterministic)),
                                                      lambda = c(new_offspring_deterministic)),
                                         nrow = nrow(new_offspring_deterministic))
      new_offspring <- apply(new_offspring_stochastic, 2, sum)
      
      population <- t(apply(survival_stochastic[seq_len(ncol(population)), , ], c(1, 2), sum))
      population[ , 1] <- population[ , 1] + new_offspring
      
    } else {
      
      population <- t(sapply(seq_len(n_cells),
                             function(x) transition_array[ , , x] %*% matrix(population[x, ])))
      
    }
    
    
    # put back in the raster
    population_raster[cell_idx] <- population
    
    landscape$population <- population_raster
    
    landscape
  }
  
  as.population_growth(pop_dynamics)
  
}


# spatial_transition <- function (survival_raster,
#                                 fecundity_raster,
#                                 demographic_stochasticity = TRUE,
#                                 global_stochasticity = 0,
#                                 local_stochasticity = 0) {
#   
#   # as above
#   
# }


# functional_transition <- function (survival_function,
#                                    fecundity_function,
#                                    demographic_stochasticity = TRUE,
#                                    global_stochasticity = 0,
#                                    local_stochasticity = 0) {
#   
#   # as above
#   
# }

# # e.g.:
# survival_fun <- function (landscape) {
#   cov <- landscape$time_since_fire
#   if (is.null(suit)) stop ()
#   
#   cell_idx <- which(!is.na(raster::getValues(cov)))
#   cov_vec <- raster::extract(cov, cell_idx)
#   
#   multipliers <- c(0.5, 0.3, 0.2)
#   cov_mat <- outer(log(cov_vec), multipliers, FUN = "*")
#   plogis(cov_mat)
# 
# }


# fecundity_fun <- function (landscape) {
#   cov <- landscape$time_since_fire
#   if (is.null(suit)) stop ()
#   
#   cell_idx <- which(!is.na(raster::getValues(cov)))
#   cov_vec <- raster::extract(cov, cell_idx)
#   
#   fec_vec <- exp(log(cov_vec) * 0.8)
#   cbind(0, 0, fec_vec)
# }



#' @rdname population_dynamics_functions
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the fast kernel-based dispersal function to modify the  
#' # population using a user-defined diffusion distribution and
#' # a fast-fourier transformation (FFT) computational algorithm:
#'
#' test_kern_dispersal <- fast_dispersal()
fast_dispersal <- function(
  dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.1),
  stages = NULL
) {
  
  pop_dynamics <- function(landscape, timestep) {
    
    # Which stages can disperse
    which_stages_disperse <- if (is.null(stages)) {
      seq_len(raster::nlayers(landscape$population))
    } else {
      stages
    }
    
    # Apply dispersal to the population
    # (need to run this separately for each stage)
    for (stage in which_stages_disperse) {
      landscape$population[[stage]][] <- dispersalFFT(
        popmat = raster::as.matrix(
          landscape$population[[stage]]
        ),
        fs = setupFFT(
          x = seq_len(raster::ncol(landscape$population)),
          y = seq_len(raster::nrow(landscape$population)),
          f = function(d) {
            disp <- dispersal_kernel(d)
            disp / sum(disp)
          }
        )
      )
    }
    
    landscape
    
  }
  
  as.population_fast_dispersal(pop_dynamics)
  
}


#' @rdname population_dynamics_functions
#' 
#' @export
#' 
#' @examples
#' 
#' # Use the probabilistic kernel-based dispersal function to modify the  
#' # population using a user-defined diffusion distribution
#' # and an arrival probability layers (e.g. habitat suitability):
#'
#' test_kern_dispersal <- kernel_dispersal()
kernel_dispersal <- function(
  distance_function = function(from, to) sqrt(rowSums(sweep(to, 2, from)^2)),
  dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.1),
  arrival_probability = c("both", "suitability", "carrying_capacity"),
  stages = NULL,
  demographic_stochasticity = FALSE
) {
  
  arrival_probability <- match.arg(arrival_probability)
  
  pop_dynamics <- function(landscape, timestep) {
    
    # check the required landscape rasters are available
    layers <- arrival_probability
    if (layers == "both")
      layers <- c("suitability", "carrying_capacity")
    
    missing_layers <- vapply(landscape[layers], is.null, FUN.VALUE = FALSE)
    if (any(missing_layers)) {
      
      missing_text <- paste(paste("a", layers[missing_layers], "raster"),
                            collapse = " and ")
      
      stop ("kernel_dispersal requires landscape to have ", missing_text,
            call. = FALSE)
    }
    
    # Which stages can disperse
    which_stages_disperse <- if (is.null(stages)) {
      seq_len(raster::nlayers(landscape$population))
    } else {
      stages      
    }
    
    # Which stages contribute to density dependence.
    which_stages_density <- if (is.null(stages)) {
      seq_len(raster::nlayers(landscape$population))
    } else {
      stages
    }
    
    # Extract locations as x and y coordinates from landscape (ncells x 2)
    xy <- raster::xyFromCell(
      landscape$population,
      seq(raster::ncell(landscape$population))
    )
    
    # Rescale locations
    xy <- sweep(xy, 2, raster::res(landscape$population), "/")
    
    delayedAssign(
      "habitat_suitability_values",
      raster::getValues(landscape$suitability)
    )
    
    delayedAssign(
      "carrying_capacity_proportion",
      raster::getValues(
        raster::calc(
          raster::stack(landscape$population)[[
            which_stages_density
            ]],
          sum
        ) /
          landscape$carrying_capacity
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
        landscape$population[[stage]]
      )
      
      # Only non-zero arrival prob cells can receive individuals
      can_arriv <- which(arrival_prob_values > 0 & !is.na(arrival_prob_values))
      
      for (stage in which_stages_disperse) {
        
        # Extract the population values
        population_values <- raster::getValues(
          landscape$population[[stage]]
        )
        
        # does cell have dispersing individuals?
        # multiply by proportion that disperses...
        # has_disperse <- which(proportion_disperses > 0 & !is.na(proportion_disperses))
        
        # Only non-zero population cells can contribute - needs to be a binomial
        # realisation of a proportion disperses function
        has_pop <- which(population_values > 0 & !is.na(population_values))
        
        contribute <- function(i) {
          # distance btw ith cell and all the cells that it can contribute to.
          contribution <- distance_function(xy[i, ], xy[can_arriv, ])
          contribution <- dispersal_kernel(contribution)
          contribution <- contribution * arrival_prob_values[can_arriv]
          # Standardise contributions and round them if demographic_stochasticity = TRUE
          contribution <- contribution / sum(contribution)
          contribution <- contribution * population_values[i]
          if (identical(demographic_stochasticity, FALSE)) return(contribution)
          contribution_int <- floor(contribution)
          idx <- utils::tail(
            order(contribution - contribution_int),
            round(sum(contribution)) - sum(contribution_int)
          )
          contribution_int[idx] <- contribution_int[idx] + 1
          contribution_int
        }
        
        landscape$population[[stage]][can_arriv] <- rowSums(
          vapply(has_pop, contribute, as.numeric(can_arriv))
        )
      }
    }
    
    landscape
    
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
                                         arrival_probability = "suitability",
                                         carrying_capacity = "carrying_capacity") {
  
  pop_dynamics <- function (landscape, timestep) {
    
    population_raster <- landscape$population
    arrival_prob <- landscape[[arrival_probability]]
    carrying_capacity <- landscape[[carrying_capacity]]
    
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
    
    landscape$population <- population_raster
    
    landscape
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
#' test_ca_dispersal <- translocation(source_layer = pop_source,
#'                                        sink_layer = pop_sink,
#'                                        stages = NULL,
#'                                        effect_timesteps = 1)

translocation <- function (source_layer, sink_layer, stages = NULL, effect_timesteps = NULL) {
  
  pop_dynamics <- function (landscape, timestep) {
    
    if (timestep %in% effect_timesteps) {
      
      population_raster <- landscape$population
      nstages <- raster::nlayers(population_raster)
      
      # get population as a matrix
      idx <- which(!is.na(raster::getValues(population_raster[[1]])))
      population_matrix <- raster::extract(population_raster, idx)
      
      source <- raster::extract(source_layer, idx)
      sink <- raster::extract(sink_layer, idx)
      
      if (is.null(stages)) stages <- seq_len(nstages)
      
      for (stage in stages) {
        
        if (any(population_matrix[ , stage] - source < 0)) {
          warning("The proposed number of translocated individuals do not exist for\nlife-stage ", stage," in timestep ", timestep, " - only the maximum number of available\nindividuals in each cell will be translocated. Please check the\nspecified source and sink layers.")
        }
        
        population_matrix[ , stage] <- population_matrix[ , stage] + sink - pmin(source, population_matrix[ , stage])
        
      }
      
      # put back in the raster
      population_raster[idx] <- population_matrix
      
      landscape$population <- population_raster
      
    }
    
    landscape
    
  }
  
  as.population_translocation(pop_dynamics)
  
}


#' @rdname population_dynamics_functions
#'
#' @export
#' 
#' @examples
#' 
#' test_pop_dd <- population_cap()

population_cap <- function (stages = NULL) {
  
  pop_dynamics <- function (landscape, timestep) {
    
    population_raster <- landscape$population
    carrying_capacity <- landscape$carrying_capacity
    
    # get population as a matrix
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population_matrix <- raster::extract(population_raster, idx)
    carrying_capacity <- raster::extract(carrying_capacity, idx)
    
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
    
    landscape$population <- population_raster 
    
    landscape
  }
  
  as.population_density_dependence(pop_dynamics)
  
}

##########################
### internal functions ###
##########################

as.population_growth <- function (simple_growth) {
  as_class(simple_growth, "population_dynamics", "function")
}

as.population_fast_dispersal <- function (fast_dispersal) {
  as_class(fast_dispersal, "population_dynamics", "function")
}

as.population_kernel_dispersal <- function (kernel_dispersal) {
  as_class(kernel_dispersal, "population_dynamics", "function")
}

as.population_ca_dispersal <- function (ca_dispersal) {
  as_class(ca_dispersal, "population_dynamics", "function")
}

as.population_translocation <- function (translocation) {
  as_class(translocation, "population_dynamics", "function")
}

as.population_density_dependence <- function (density_dependence) {
  as_class(density_dependence, "population_dynamics", "function")
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