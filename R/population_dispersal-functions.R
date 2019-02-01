#' @useDynLib steps
#' @importFrom Rcpp sourceCpp
NULL
#' How the population disperses in a landscape.
#'
#' Pre-defined functions to define population dispersal during a simulation.
#'
#' @name population_dispersal_functions
#'
#' @param demographic_stochasticity should demographic stochasticity be used in
#'   dispersal?
#' @param dispersal_kernel a single or list of user-defined distance dispersal
#'   kernel functions
#' @param dispersal_proportion proportions of individuals (0 to 1) that can
#'   disperse in each life stage
#' @param arrival_probability a raster layer that controls where individuals can
#'   disperse to (e.g. habitat suitability)
#' @param dispersal_distance the maximum distance (in cells) that each life stage
#'   can disperse (in kernel-based dispersal this truncates the dispersal curve)
#' @param barrier_type if barrier map is used, does it stop (0 - default) or
#'   kill (1) individuals
#' @param dispersal_steps number of dispersal steps to take before stopping
#' @param use_barriers should dispersal barriers be used? If so, a barriers map
#'   must be provided
#' @param barriers_map a raster layer that contains cell values of 0 (no
#'   barrier) and 1 (barrier)
#' @param carrying_capacity a raster layer that specifies the carrying capacity
#'   in each cell
#'
#' @rdname population_dispersal_functions
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
  dispersal_proportion = 1,
  demographic_stochasticity = FALSE
) {
  
  pop_dynamics <- function(landscape, timestep) {
    
    n_stages <- raster::nlayers(landscape$population)
    
    #poptot <- sum(raster::cellStats(landscape$population, sum))

    if (length(dispersal_proportion) < n_stages) {
      if (timestep == 1) cat ("    ", n_stages, "life stages exist but", length(dispersal_proportion),"dispersal proportion(s) of", dispersal_proportion,"were specified. Is this what was intended?")
      dispersal_proportion <- rep(dispersal_proportion, n_stages)[1:n_stages]
    }
    
    # Which stages can disperse
    which_stages_disperse <- which(dispersal_proportion > 0)
    
    # Apply dispersal to the population
    # (need to run this separately for each stage)
    for (stage in which_stages_disperse) {
      
      pop <- landscape$population[[stage]]
      
      pop_dispersing <- landscape$population[[stage]] * dispersal_proportion[[stage]]
      pop_staying <- pop - pop_dispersing
      
      pop_dispersed <- dispersalFFT(
        popmat = raster::as.matrix(
          pop_dispersing
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
      
      raster::values(pop_dispersing) <- pop_dispersed
      pop <- pop_staying + pop_dispersing
      landscape$population[[stage]] <- pop
      
    }
    
    #cat("Pre-Post Population:", poptot, sum(raster::cellStats(landscape$population, sum)), "(Timestep", timestep, ")", "\n")
    
    landscape
    
  }
  
  as.population_fast_dispersal(pop_dynamics)
  
}


#' @rdname population_dispersal_functions
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
  dispersal_kernel = exponential_dispersal_kernel(distance_decay = 10000),
  dispersal_distance = Inf,
  arrival_probability = c("both", "suitability", "carrying_capacity"),
  dispersal_proportion = 1,
  demographic_stochasticity = TRUE
) {
  
  arrival_probability <- match.arg(arrival_probability)
  
  pop_dynamics <- function(landscape, timestep) {
    
    # extract locations as x and y coordinates from landscape (ncells x 2)
    xy <- raster::xyFromCell(
      landscape$population,
      seq(raster::ncell(landscape$population))
    )
    
    # rescale locations
    xy <- sweep(xy, 2, raster::res(landscape$population), "/")

    distance_info <- steps_stash$distance_info
    
    # if (is.null(dist_list)) {
    #   dist_mat <- as.matrix(dist(xy))
    #   dist_mat[dist_mat > dispersal_distance] <- NA
    #   dist_list <- apply(dist_mat, 1, function(x) which(!is.na(x)))
    # }
    
    if (is.null(distance_info)) {
      distance_info <- get_distance_info(res = raster::res(landscape$population),
                                                      max_distance = dispersal_distance)
    }

    # how many life stages?
    n_stages <- raster::nlayers(landscape$population)

    # check the required landscape rasters/functions are available
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

    # get non-NA cells
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
    
    if (exists("carrying_capacity_function", envir = steps_stash)) {
      if (raster::nlayers(landscape$suitability) > 1) landscape$carrying_capacity[cell_idx] <- steps_stash$carrying_capacity_function(landscape$suitability[[timestep]][cell_idx])
      else landscape$carrying_capacity[cell_idx] <- steps_stash$carrying_capacity_function(landscape$suitability[cell_idx])
    }
    
    if (length(dispersal_proportion) < n_stages) {
      if (timestep == 1) cat ("    ", n_stages, "life stages exist but", length(dispersal_proportion),"dispersal proportion(s) of", dispersal_proportion,"were specified. Is this what was intended?")
      dispersal_proportion <- rep(dispersal_proportion, n_stages)[1:n_stages]
    }
    
    # Which stages can disperse
    which_stages_disperse <- which(dispersal_proportion > 0)
    
    # Which stages contribute to density dependence.
    which_stages_density <- which(dispersal_proportion > 0)
    
    delayedAssign(
      "habitat_suitability_values",
      if (raster::nlayers(landscape$suitability) > 1) raster::getValues(landscape$suitability[[timestep]])
      else raster::getValues(landscape$suitability)
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
      both = habitat_suitability_values * (1 - carrying_capacity_proportion),
      suitability = habitat_suitability_values,
      carrying_capacity = (1 - carrying_capacity_proportion)
    )
    
    # Only non-zero arrival prob cells can receive individuals
    can_arriv <- which(arrival_prob_values > 0 & !is.na(arrival_prob_values))

    for (stage in which_stages_disperse) {
      
      # Extract the population values
      pop <- landscape$population[[stage]]
      
      # Calculate the proportion dispersing
      pop_dispersing <- raster::getValues(
        landscape$population[[stage]] * dispersal_proportion[[stage]]
      )
      
      # Calculate the proportion not dispersing
      pop_staying <- pop - pop_dispersing

      # does cell have dispersing individuals?
      # multiply by proportion that disperses...
      # has_disperse <- which(proportion_disperses > 0 & !is.na(proportion_disperses))
      
      # Only non-zero population cells can contribute - needs to be a binomial
      # realisation of a proportion disperses function
      has_pop <- which(pop_dispersing > 0 & !is.na(pop_dispersing))

      contribute <- function(i) {
        
        # distance btw ith cell and all the cells that it can contribute to.
        #dists <- sqrt(rowSums(sweep(xy[can_arriv, ], 2, xy[i, ])^2))
        
        # which cells are within the maximum distance
        #in_range <- which(dists <= dispersal_distance)

        # what are the probabilities of dispersing given the distances
        #contribution <- rep(0, length(dists))
        #contribution[in_range] <- dispersal_kernel(dists[in_range])
        
        # what cells can be dispersed to based on the distance cap
        # within_dist <- which(round(local_dist_matrix[[i]][[1]], 3) <= dispersal_distance)
        
        # arriv_dist <- intersect(within_dist, can_arriv)
        
        #contribution <- distance_function(xy[i, ], xy[can_arriv, ])
        
        destinations <- get_ids_dists(cell_id = i,
                               distance_info = distance_info,
                               raster_dim = dim(landscape$population[[1]])[-3])
        
        destination_ids <- destinations[, 1]
        destination_dists <- destinations[, 2]

        # account for arrival posibility
        keep_destination <- destination_ids %in% intersect(can_arriv, destination_ids)
        destination_ids <- destination_ids[keep_destination]
        destination_dists <- destination_dists[keep_destination]
        
        contribution <- dispersal_kernel(destination_dists)
        
        # probability of dispersing multiplied by probability of arrival 
        contribution <- contribution * arrival_prob_values[destination_ids]
        
        # standardise contributions
        contribution <- contribution / sum(contribution)
        
        # estimate population dispersing based on contribution
        contribution <- contribution * pop_dispersing[i]
        
        final_pop <- rep(0, length(can_arriv))
        final_pop[destination_ids] <- contribution
        
        # round contributions if demographic_stochasticity = TRUE
        if (identical(demographic_stochasticity, FALSE)) return(final_pop)

        final_pop_int <- floor(final_pop)
        idx <- utils::tail(
          order(final_pop - final_pop_int),
          round(sum(final_pop)) - sum(final_pop_int)
        )
        final_pop_int[idx] <- final_pop_int[idx] + 1
        final_pop_int
        ####################################
      }
      
      landscape$population[[stage]][can_arriv] <- rowSums(
        vapply(has_pop, contribute, as.numeric(can_arriv))
      ) + pop_staying[can_arriv]
      
      #values(pop_dispersing) <- pop_dispersed
      #pop <- pop_staying + pop_dispersing
      #landscape$population[[stage]] <- pop
      
      
    }
    
    #cat("Pre-Post Population:", poptot, sum(raster::cellStats(landscape$population, sum)), "(Timestep", timestep, ")", "\n")
    
    landscape
    
  }
  
  as.population_kernel_dispersal(pop_dynamics)
  
}


#' @rdname population_dispersal_functions
#' 
#' @export
#' 
#' @examples
#'
#' # Use the cellular automata dispersal function to modify  
#' # the population using rule-based cell movements:
#' 
#' test_ca_dispersal <- cellular_automata_dispersal()

cellular_automata_dispersal <- function (dispersal_distance = 1,
                                         dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.1),
                                         dispersal_proportion = 1,
                                         barrier_type = 0,
                                         dispersal_steps = 1,
                                         use_barriers = FALSE,
                                         barriers_map = NULL,
                                         arrival_probability = "suitability",
                                         carrying_capacity = "carrying_capacity") {
  
  pop_dynamics <- function (landscape, timestep) {
    
    population_raster <- landscape$population
    
    # Get non-NA cells
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    
    if (exists("carrying_capacity_function", envir = steps_stash)) {
      if (raster::nlayers(landscape$suitability) > 1) landscape$carrying_capacity[idx] <- steps_stash$carrying_capacity_function(landscape$suitability[[timestep]][idx])
      else landscape$carrying_capacity[idx] <- steps_stash$carrying_capacity_function(landscape$suitability[idx])
    }
    
    if (raster::nlayers(landscape$suitability) > 1 & arrival_probability == "suitability") arrival_prob <- landscape[[arrival_probability]][[timestep]]
    else arrival_prob <- landscape[[arrival_probability]]
    
    carrying_capacity <- landscape[[carrying_capacity]]
    
    #poptot <- sum(raster::cellStats(landscape$population, sum))
    
    # get population as a matrix
    population <- raster::extract(population_raster, idx)

    n_stages <- raster::nlayers(population_raster)
    
    if (length(dispersal_proportion) < n_stages) {
      if (timestep == 1) cat ("    ", n_stages, "life stages exist but", length(dispersal_proportion),"dispersal proportion(s) of", dispersal_proportion,"were specified. Is this what was intended?")
      dispersal_proportion <- rep(dispersal_proportion, n_stages)[1:n_stages]
    }

    if (length(dispersal_distance) < n_stages) {
      if (timestep == 1) cat ("    ", n_stages, "life stages exist but", length(dispersal_distance),"dispersal distance(s) of", dispersal_distance,"were specified. Is this what was intended?")
      dispersal_distance <- rep(dispersal_distance, n_stages)[1:n_stages]
    }
    
    dispersal_vector <- lapply(dispersal_distance, function(x) dispersal_kernel(seq_len(x)-1))
        
    # identify dispersing stages
    which_stages_disperse <- which(dispersal_proportion > 0)
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
    
    # could do this in parallel
    for (i in which_stages_disperse){
      population_raster[[i]][] <- rcpp_dispersal(raster::as.matrix(population_raster[[i]]),
                                                 raster::as.matrix(carrying_capacity),
                                                 raster::as.matrix(arrival_prob),
                                                 raster::as.matrix(barriers_map),
                                                 as.integer(barrier_type),
                                                 use_barriers,
                                                 as.integer(dispersal_steps),
                                                 as.integer(dispersal_distance[i]),
                                                 as.numeric(unlist(dispersal_vector[i])),
                                                 as.numeric(dispersal_proportion[i])
      )$dispersed_population
      
    }
    
    landscape$population <- population_raster
    
    #cat("Pre-Post Population:", poptot, sum(raster::cellStats(landscape$population, sum)), "(Timestep", timestep, ")", "\n")
    
    landscape
  }
  
  as.population_ca_dispersal(pop_dynamics)
  
}


##########################
### internal functions ###
##########################

as.population_fast_dispersal <- function (fast_dispersal) {
  as_class(fast_dispersal, "population_dynamics", "function")
}

as.population_kernel_dispersal <- function (kernel_dispersal) {
  as_class(kernel_dispersal, "population_dynamics", "function")
}

as.population_ca_dispersal <- function (ca_dispersal) {
  as_class(ca_dispersal, "population_dynamics", "function")
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
  
  # duplicate popmat to create 'before' condition
  popmat_orig <- popmat
  
  # check for missing values and replace with zeros
  if (any(is.na(popmat))) {
    popmat[is.na(popmat)] <- 0
  }
    
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
  
  # get proportion of population that dispersed into NA areas
  prop_out <- sum(pop_new[is.na(popmat_orig)]) / sum(pop_new[!is.na(popmat_orig)])

  # check for NaN and replace with zero
  if (is.nan(prop_out)) prop_out <- 0
    
  # return NA values to matrix
  pop_new[is.na(popmat_orig)] <- NA
  
  # increase all non-NA cells by inverse of proportion in NA areas
  pop_new[!is.na(popmat_orig)] <- pop_new[!is.na(popmat_orig)] * (1 + prop_out)

  # make sure none are lost or gained (unless all are zeros)
  if (any(pop_new[!is.na(popmat_orig)] > 0)) {
    pop_new[!is.na(popmat_orig)] <- stats::rmultinom(1, size = sum(popmat), prob = pop_new[!is.na(popmat_orig)])    
  }
  
  pop_new
}

seq_range <- function (range, by = 1) seq(range[1], range[2], by = by)

# compute the *relative* distances to all neightbouring cells within a maximum
# distance
get_distance_info <- function(res, max_distance) {
  
  width <- ceiling(max_distance / min(res)) + 1
  id <- (1:width) - 1
  xy <- expand.grid(id * res[1], id * res[2])
  cell_coord <- expand.grid(id, id)
  cell_coord <- as.matrix(cell_coord)
  colnames(cell_coord) <- NULL
  dists <- as.matrix(dist(xy))[1, ]
  keep <- dists < max_distance
  
  # relative coordinates of cells that are within the distance
  ur <- cell_coord[keep, ]
  ul <- cbind(-ur[, 1], ur[, 2])
  ll <- -ur
  lr <- cbind(ur[, 1], -ur[, 2])
  
  # all coordinates in the circle
  coords <- rbind(ur, ul, ll, lr)
  
  # add on their distances
  cbind(coords, rep(dists[keep], 4))
  
}

# given a cell id, find the coordinates (in number of cells from the origin)
id2coord <- function (id, dim) {
  # get coordinates in cell numbers. 
  # index from 0
  id0 <- id - 1
  
  # how many rows have been covered
  row0 <- id0 %/% dim[1]
  # how many columns have been covered in the most recent row
  col0 <- id0 %% dim[1]
  # combine and switch back to indexing from 1
  coord <- cbind(col0, row0) + 1
  # make sure they are within the raster (for completeness)
  max_id <- prod(dim)
  invalid_id <- id > max_id | id < 1
  coord[invalid_id, ] <- NA
  coord
  
}

# returns NAs where the coords are outside the raster. coord must be a 2-column
# matrix
coord2id <- function (coord, dim) {
  start_of_row <- (coord[, 2] - 1) * dim[2]
  id <- start_of_row + coord[, 1]
  
  max_id <- prod(dim)
  invalid_id <- id > max_id | id < 1
  id[invalid_id] <- NA
  id
  
}

# given a cell number, compute the cell ids that are within the maximum
# distance, and get the distances to those
get_ids_dists <- function(cell_id, distance_info, raster_dim) {
  
  # get origin coordinate  
  origin_coord <- id2coord(cell_id, raster_dim)  
  
  # relative coordinates (in numbers of cells) of cells within max distance, and
  # their distances from this cell
  rel_coords <- distance_info[, 1:2]
  dists <- distance_info[, 3]
  
  # compute absolute coordinates
  abs_coords <- sweep(rel_coords, 2, origin_coord, "+") 
  
  cell_ids <- coord2id(abs_coords, raster_dim)
  valid <- !is.na(cell_ids)
  
  # return a 2-column matrix of cells id and distances for all cells within the
  # maximum distances, and within the raster
  cbind(cell_ids[valid], dists[valid])
  
}

