#' @useDynLib steps
#' @importFrom Rcpp sourceCpp
NULL

#' How the population disperses in a landscape.
#'
#' Pre-defined or custom functions to define population dispersal during a
#' simulation. Each dispersal method uses different computing resources
#' and may be applicable to different simulation scenarios.
#' 
#' The fast_dispersal function uses kernel-based dispersal
#' to modify the population with a user-defined diffusion distribution
#' and a fast-fourier transformation (FFT) computational algorithm. It
#' is computationally efficient and very fast, however, only useful for
#' situations where dispersal barriers or arrival based on habitat or
#' carrying capacity are not required. In other words, organisms
#' can disperse in all directions and to all cells in the landscape.
#' 
#' The kernel_dispersal function employs a probabilistic
#' kernel-based dispersal algorithm to modify the population
#' using a user-defined diffusion distribution, arrival
#' probability layers (e.g. habitat suitability), and growth
#' limiting layers (e.g. carrying capacity). This function is much
#' slower than the fast_dispersal, however, respects dispersal
#' limitations which may be more ecologically appropriate. Further, 
#' the kernel-based dispersal function utilises a mechanism to
#' optimise computational performance in which it switches between
#' pre-allocating cell movements based on the available memory of
#' the host computer (faster but more memory intensive) or executing
#' cell movements in sequence (slower but less memory intensive).
#' 
#' The cellular_automata_dispersal function modifies populations
#' using rule-based cell movements. This function allows the use
#' of barriers in the landscape to influence dispersal. The
#' function is computationally efficient, however, because
#' individuals are dispersed, performance scales with the
#' population sizes in each cell across a landscape.
#'
#' @name population_dispersal_functions
#'
#' @param dispersal_kernel a single built-in or user-defined distance dispersal
#'   kernel function.
#' @param dispersal_proportion proportions of individuals (0 to 1) that can
#'   disperse in each life stage. Can also be a custom function that relates
#'   the proportion dispersing to a spatial layer in the landscape object
#'   (e.g. carrying capacity).
#' @param arrival_probability the name of a spatial layer in the landscape object
#'   that controls where individuals can disperse to (e.g. "suitability") or
#'   "none" to allow individuals to disperse to all non-NA cells.
#' @param carrying_capacity the name of a spatial layer in the landscape object
#' that specifies the carrying capacity in each cell.
#' @param max_distance the maximum distance that each life stage can
#'   disperse in spatial units of the landscape (in kernel-based dispersal
#'   this truncates the dispersal curve).
#' @param max_cells the maximum number of cell movements that each individual in
#'   each life stage can disperse in whole integers.
#' @param barriers_map the name of a spatial layer in the landscape object that
#'   contains cell values between 0 (no barrier) and 1 (full barrier) Any
#'   values between 0 and 1 indicate the permeability of the barrier.   
#' @param use_suitability should habitat suitability be used to control the
#'   likelihood of individuals dispersing into cells? The default is TRUE. Note,
#'   if a barrier map is also provided, the suitability map is multiplied with
#'   the barrier map to generate a permeability map of the landscape.
NULL

#' @rdname population_dispersal_functions
#' 
#' @export
#' 
#' @examples
#' 
#' test_fast_dispersal <- fast_dispersal()

fast_dispersal <- function(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 0.1),
                           dispersal_proportion = all_dispersing()) {
  
  pop_dynamics <- function(landscape, timestep) {
    
    n_stages <- raster::nlayers(landscape$population)
    
    dispersal_proportion <- dispersal_proportion(landscape, timestep)
    
    # Which stages can disperse
    which_stages_disperse <- which(dispersal_proportion > 0)
    
    # Apply dispersal to the population
    # (need to run this separately for each stage)
    for (stage in which_stages_disperse) {
      
      pop <- landscape$population[[stage]]
      
      pop_dispersing <- landscape$population[[stage]] * dispersal_proportion[[stage]]
      pop_staying <- pop - pop_dispersing
      
      # round population staying
      idx <- not_missing(pop_staying)
      pop_staying_vec <- raster::extract(pop_staying, idx)
      pop_staying_vec <- round_pop(pop_staying_vec)
      pop_staying[idx] <- pop_staying_vec
      
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
      
      pop_dispersing[] <- pop_dispersed
      pop <- pop_staying + pop_dispersing
      landscape$population[[stage]] <- pop
      
    }
    
    #cat("Pre-Post Population:", poptot, sum(raster::cellStats(landscape$population, sum)), "(Timestep", timestep, ")", "\n")
    
    landscape
    
  }
  
  as.population_dispersal(pop_dynamics)
  
}


#' @rdname population_dispersal_functions
#' 
#' @export
#' 
#' @importFrom memuse Sys.meminfo mu.size
#' 
#' @examples
#' 
#' test_kern_dispersal <- kernel_dispersal()

kernel_dispersal <- function (dispersal_kernel = exponential_dispersal_kernel(distance_decay = 1),
                              max_distance = Inf,
                              arrival_probability = c("both", "suitability", "carrying_capacity", "none"),
                              dispersal_proportion = all_dispersing()) {
  
  arrival_probability <- match.arg(arrival_probability)
  
  pop_dynamics <- function(landscape, timestep) {
    
    # if max_distance is the default (Infinite), rescale to maximum distance of landscape
    default_distance <- identical(max_distance, Inf)
    
    if (default_distance) {
      n_rows <- raster::nrow(landscape$population[[1]])
      n_cols <- raster::ncol(landscape$population[[1]])
      res <- raster::res(landscape$population[[1]])
      max_distance <- sqrt( (n_cols * res[1])^2 + (n_rows * res[2])^2 )
    }
    
    
    distance_list <- steps_stash$distance_list
    if (is.null(distance_list)) {
      # what are dimensions of raster (lazyeval was causing this to be rerun on
      # every contribute() call! so force execution of this here)
      raster_dim <- dim(landscape$population[[1]])[-3]
      raster_dim <- force(raster_dim)
      
      distance_info <- get_distance_info(res = raster::res(landscape$population),
                                         max_distance = max_distance)
      
      sys_mem_available <- memuse::Sys.meminfo()$freeram
      
      sys_mem_available <- memuse::mu.size(sys_mem_available, as.is = FALSE) * 0.8
      
      n_elem <- nrow(distance_info) * raster::ncell(landscape$population)
      sys_mem_required <- (n_elem * (64 + 32)) / 8
      
      if (sys_mem_required < sys_mem_available) {
        distance_list <- steps_stash$distance_list <- lapply(seq_len(raster::ncell(landscape$population)),
                                                             function (x) get_ids_dists(cell_id = x,
                                                                                        distance_info = distance_info,
                                                                                        raster_dim = raster_dim))
      }
      
    }
    
    # how many life stages?
    n_stages <- raster::nlayers(landscape$population)
    
    # create masking layer
    mask <- raster::getValues(landscape$population[[1]])
    mask[!is.na(mask)] <- 1
    
    # check the required landscape rasters/functions are available
    layers <- arrival_probability
    if (layers == "both") {
      layers <- c("suitability", "carrying_capacity")
      
      missing_layers <- vapply(landscape[layers], is.null, FUN.VALUE = FALSE)
      if (any(missing_layers)) {
        
        missing_text <- paste(paste("a", layers[missing_layers], "raster"),
                              collapse = " and ")
        
        stop ("kernel_dispersal requires landscape to have ", missing_text,
              call. = FALSE)
      }
    }
    
    # find out which stages contribute to density. population dynamics will have
    # copied this information over from the population density dependence
    # function, if it was specified. Otherwise, use all stages
    if (!exists("density_dependence_stages")) {
      density_dependence_stages <- seq_len(raster::nlayers(landscape$population))
    }
    
    dispersal_proportion <- dispersal_proportion(landscape, timestep)
    
    # Which stages can disperse
    which_stages_disperse <- which(dispersal_proportion > 0)

    # get non-NA cells
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
    
    delayedAssign(
      "habitat_suitability_values",
      if (raster::nlayers(landscape$suitability) > 1) raster::getValues(landscape$suitability[[timestep]])
      else raster::getValues(landscape$suitability)
    )
    
    if ("carrying_capacity" %in% layers) {
      
      cc <- get_carrying_capacity(landscape, timestep)
      
      if (is.null(cc)) {
        stop ("carrying capacity must be specified in the landscape object to use carrying capacity arrival probabilities",
              call. = FALSE)
      }

      delayedAssign(
        "carrying_capacity_proportion",
        raster::getValues(
          raster::calc(
            raster::stack(landscape$population)[[
              density_dependence_stages
              ]],
            sum
          ) / cc
        )
      )
      
      cc_values <- raster::getValues(cc)
      
    } else {
      
      cc_values <- NULL
      density_dependence_stages <- NULL
      
    }
    
    arrival_prob_values <- switch(
      arrival_probability,
      both = habitat_suitability_values * (1 - carrying_capacity_proportion),
      suitability = habitat_suitability_values,
      carrying_capacity = (1 - carrying_capacity_proportion),
      none = mask
    )
    
    # Only non-zero arrival prob cells can receive individuals
    can_arriv_ids <- which(arrival_prob_values > 0 & !is.na(arrival_prob_values))
    
    # Extract the population values
    pop <- raster::getValues(landscape$population)
    
    # Only non-zero population cells can contribute
    #pop_dispersing_vec <- pop_dispersing[]
    has_pop_ids <- which(rowSums(pop) > 0 & !is.na(pop[, 1]))
    
    # store the original population, so that individuals don't disperse twice
    original_pop <- pop

    for (origin in sample(has_pop_ids)) {
      for (stage in sample(which_stages_disperse)) {
        pop <- disperse(origin = origin,
                        stage = stage,
                        pop = pop,
                        original_pop = original_pop,
                        prop_dispersing = dispersal_proportion,
                        can_arriv_ids = can_arriv_ids,
                        arrival_prob_values = arrival_prob_values,
                        dispersal_kernel = dispersal_kernel,
                        carrying_capacity = cc_values,
                        density_dependence_stages = density_dependence_stages,
                        distance_list = distance_list,
                        distance_info = distance_info,
                        raster_dim = raster_dim)
      }
    }

    landscape$population[] <- pop
    
    #cat("Pre-Post Population:", poptot, sum(raster::cellStats(landscape$population, sum)), "(Timestep", timestep, ")", "\n")
    
    landscape
    
  }
  
  as.population_dispersal(pop_dynamics)
  
}

#' @rdname population_dispersal_functions
#'
#' @export
#' 
#' @examples
#'
#' test_ca_dispersal <- cellular_automata_dispersal()

cellular_automata_dispersal <- function (max_cells = Inf,
                                         dispersal_proportion = all_dispersing(),
                                         barriers_map = NULL,
                                         use_suitability = TRUE,
                                         carrying_capacity = "carrying_capacity") {
  
  # are there suitability and carrying_capacity landscape objects specified?
  carrying_cap <- identical(carrying_capacity, "carrying_capacity")
  
  pop_dynamics <- function (landscape, timestep) {
    
    # read population raster stack from landscape
    population_raster <- landscape$population
    
    # get non-NA cells
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    
    # is suitability specified without raster existing in landscape?
    if (use_suitability & is.null(landscape$suitability)) {
      stop("A habitat suitability raster is missing from the landscape object")
    }
    
    # handle input suitability raster stacks
    if (use_suitability & raster::nlayers(landscape$suitability) > 1) {
      suitability_map <- landscape$suitability[[timestep]]
    } else {
      suitability_map <- landscape$suitability
    }
    
    # is carrying_capacity specified without raster existing in landscape?
    if (carrying_cap & is.null(landscape$carrying_capacity)) {
      stop("A carrying capacity object is missing from the landscape object")
    }
    
    # handle carrying_capacity as raster, function, or other spatial object in
    # landscape
    if (carrying_cap) {
      cc <- get_carrying_capacity(landscape, timestep)
    } else {
      cc <- landscape[[carrying_capacity]]
    }
    
    # get population as a matrix
    population <- raster::extract(population_raster, idx)
    
    # get number of life-stages
    n_stages <- raster::nlayers(population_raster)
    
    # work out dispersal proportions for eligible life-stages with input function
    dispersal_proportion <- dispersal_proportion(landscape, timestep)
    
    # if max_cells is the default (Infinite), rescale to maximum distance of landscape
    default_distance <- identical(max_cells, Inf)
    
    if (default_distance) {
      n_rows <- raster::nrow(population_raster[[1]])
      n_cols <- raster::ncol(population_raster[[1]])
      res <- raster::res(population_raster[[1]])
      max_cells <- sqrt( (n_cols * res[1])^2 + (n_rows * res[2])^2 )
    }
    
    # handle dispersal distances as both scalars and vectors
    if(!default_distance) {
    warn_once(length(max_cells) < n_stages | length(max_cells) > n_stages,
              paste(n_stages,
                    "life stages exist but",
                    length(max_cells),
                    "maximum cell movement(s) of",
                    paste(max_cells, collapse = ", "),
                    "were specified.\nAll life stages will use this distance."),
              warning_name = "dispersal_distances")
    }
    
    if (length(max_cells) < n_stages | length(max_cells) > n_stages) {
      max_cells <- rep_len(max_cells, n_stages)
    }

    # identify dispersing stages
    which_stages_disperse <- which(dispersal_proportion > 0 & max_cells > 0)
    
    # if no barrier map is specified, create a barriers matrix with all zeros.
    if (is.null(barriers_map)) {
      barriers_map <- raster::calc(population_raster[[1]],
                                   function(x){x[!is.na(x)] <- 0; return(x)})
    } else {
      barriers_map <- landscape[[barriers_map]]
    }
    
    if (use_suitability) {
      permeability_map <- suitability_map * (1 - barriers_map)
    } else {
      permeability_map <- 1 - barriers_map
    }

    steps_stash$dispersal_stats_success <- replicate(n_stages, NULL, simplify = FALSE)
    steps_stash$dispersal_stats_failure <- replicate(n_stages, NULL, simplify = FALSE)
    
    # could do this in parallel
    for (i in which_stages_disperse){
      dispersed <- rcpp_dispersal(raster::as.matrix(population_raster[[i]]),
                                  raster::as.matrix(cc),
                                  raster::as.matrix(permeability_map),
                                  as.integer(max_cells[i]),
                                  as.numeric(dispersal_proportion[i]))
      
      population_raster[[i]][] <- dispersed$future_population

      steps_stash$dispersal_stats_success[[i]] <- dispersed$dispersed
      steps_stash$dispersal_stats_failure[[i]] <- dispersed$failed
    }

    landscape$population <- population_raster
    
    #cat("Pre-Post Population:", poptot, sum(raster::cellStats(landscape$population, sum)), "(Timestep", timestep, ")", "\n")
    
    landscape
  }
  
  as.population_dispersal(pop_dynamics)
  
}


##########################
### internal functions ###
##########################

as.population_dispersal <- function (dispersal) {
  as_class(dispersal, "population_dispersal", "function")
}

# as.population_fast_dispersal <- function (fast_dispersal, info = NULL) {
#   as_class(fast_dispersal, "population_fast_dispersal", "function", info = info)
# }
# 
# as.population_kernel_dispersal <- function (kernel_dispersal, info = NULL) {
#   as_class(kernel_dispersal, "population_kernel_dispersal", "function", info = info)
# }
# 
# as.population_ca_dispersal <- function (ca_dispersal, info = NULL) {
#   as_class(ca_dispersal, "population_ca_dispersal", "function", info = info)
# }

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
  
  missing <- is.na(popmat)
  # check for missing values and replace with zeros
  if (any(missing)) {
    popmat[missing] <- 0
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
  
  # return NA values to matrix
  pop_new[missing] <- NA
  
  # get proportion of population that dispersed into valid (non-NA, inside the
  # plane) cells

  prop_in <- sum(pop_new, na.rm = TRUE) / sum(popmat_orig, na.rm = TRUE)
  
  # increase all non-NA cells by inverse of proportion in valid cells
  pop_new[!missing] <- pop_new[!missing] / prop_in
  
  # make sure none are lost or gained (unless all are zeros)
  if (any(pop_new[!missing] > 0)) {
    pop_new[!missing] <- round_pop(pop_new[!missing])
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
  dists <- as.matrix(stats::dist(xy))[1, ]
  keep <- dists < max_distance
  
  # relative coordinates of cells that are within the distance
  ur <- cell_coord[keep, ]
  ul <- cbind(-ur[, 1], ur[, 2])
  ll <- -ur
  lr <- cbind(ur[, 1], -ur[, 2])
  
  # all coordinates in the circle
  coords <- rbind(ur, ul, ll, lr)
  
  # add on their distances & remove duplicates
  unique(cbind(coords, rep(dists[keep], 4)))
  
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
  
  # find coordinates outside the raster
  bad_row <- coord[, 2] < 1 | coord[, 2] > dim[1]
  bad_col <- coord[, 1] < 1 | coord[, 1] > dim[2]
  invalid <- bad_row | bad_col
  
  id[invalid] <- NA
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

disperse <- function (origin,
                      stage,
                      pop,
                      original_pop,
                      prop_dispersing,
                      can_arriv_ids,
                      arrival_prob_values,
                      dispersal_kernel,
                      carrying_capacity = NULL,
                      density_dependence_stages = NULL,
                      distance_list = NULL,
                      distance_info = NULL,
                      raster_dim = NULL) {
  
  if (is.null(distance_list)) {
    
    destinations <- get_ids_dists(cell_id = origin,
                                  distance_info = distance_info,
                                  raster_dim = raster_dim)
  } else {
    destinations <- distance_list[[origin]]
  }
  
  destination_ids <- destinations[, 1]
  destination_dists <- destinations[, 2]
  
  # index to cells in raster
  valid_arriv <- destination_ids[match(can_arriv_ids, destination_ids, 0L)]
  
  # index to cells in can_arriv_ids
  arrival_index <- match(valid_arriv, can_arriv_ids)
  
  # account for arrival possibility
  keep_destination <- destination_ids %in% valid_arriv
  destination_ids <- destination_ids[keep_destination]
  destination_dists <- destination_dists[keep_destination]
  
  prob <- dispersal_kernel(destination_dists)
  
  # probability of dispersing multiplied by probability of arrival
  prob <- prob * arrival_prob_values[destination_ids]
  
  # standardise probabilities
  prob <- prob / sum(prob)
  
  # get number dispersing and staying
  # (if this is not the first pcell considered, we use the original population
  # to make sure new arrivals don't disperse again)
  n_total <- original_pop[origin, stage]
  n_dispersing <- round_pop(n_total * prop_dispersing[stage])
  n_staying <- n_total - n_dispersing
  
  # update pop and original_pop to remove the dispersers
  new_arrivals <- pop[origin, stage] - n_total
  pop[origin, stage] <- n_staying + new_arrivals
  
  # propose some dispersals
  dispersals <- round_pop(n_dispersing * prob)

  # if we're using carrying capacity, return individuals to the origin if there's no space
  if (!is.null(carrying_capacity)) {
    
    # get space left in each, and excess dispersers
    effective_populations <- pop[destination_ids, ]
    effective_populations[, -density_dependence_stages] <- 0
    effective_population <- rowSums(effective_populations)
    space_remaining <- carrying_capacity[destination_ids] - effective_population
    
    # if there's any stochastcity, need to round down
    if (steps_stash$demo_stochasticity != "none") {
      space_remaining <- floor(space_remaining)
    }
    
    excess <- pmax(0, dispersals - space_remaining)
    destination_is_origin <- which(destination_ids == origin)
      
    # remove the excess from the dispersers, and say they are dispersing back to
    # the origin
    dispersals <- dispersals - excess
    dispersals[destination_is_origin] <- dispersals[destination_is_origin] + sum(excess)
      
  }
    
  # assign them to their population and return
  pop[destination_ids, stage] <- pop[destination_ids, stage] + dispersals
  
  pop
  
}
