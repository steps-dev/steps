#' @useDynLib dhmpr
#' @importFrom Rcpp sourceCpp
NULL

#' @title dispersal class for meta-populations
#' @rdname dispersal-class
#' @name as.dispersal
#' @description creates a function that governs dispersal capacity of life-history stages of a species. The input is a list which contains the first list is a dispersal kernel \code{dispersal} value for life-history stage, the second named list \code{} is the proportion of that stage will disperse. For example a probability of 0.2 for stage larvae means a random 20% of larve will a try and disperse to patches, the distance they can disperse is governed by the dispersal kernel (alpha). If params = NULL, a dispersal kernel of 1 is given to all stages, and all stages will attempt to undertake dispersal. If NULL is provided to the as.dispersal function, diffuse dispersal will be used based on the fast fourier transformation method ('fft').  
#' @param params A list of named lists which contain the parameters form dispersal behaviour and other parameters for dispersal modules - see details below for more information.
#' @details text describing parameter inputs in more detail.
#' \itemize{
#'  \item{"dispersal_distance"}{ list The number of cells that each life-history can disperse.}
#'  \item{"dispersal_kernel"}{ list The dispersal kernel for each life history stage. Needs to be a numeric vector that
#'  matches the length of dispersal_distance.}
#'  \item{"dispersal_proportion"}{ list The proportion of the population in each cell that will disperse. e.g 0.6 = 60\%.}
#'  \item{"dispersal_steps"} { int The number of cellular automata iterations to do in each stage. Note to self, this could be linked to time functions, ie, daily dispersal (7 dispersal steps) to match weekly fire model.}
#'  \item{"use_barriers"}{ bool To use barriers in dispersal step, only applicable to cellular automata}
#'  \item{"barrier_type"}{ type Barriers can be "blocking" or "stopping". If a barriers is blocking it will stop dispersal to that cell, but allow dispersal to other nearby cells, if they meet all the conditions of dispersal. Stopping will stop dispersal if a cell a barrier is contacted.}
#' }
#' @param method Can be either fast fourier transformation = 'fft'; or cellular automata = 'ca'. 
#' 
#' 
#' 
#' @export
#' @examples 
#' dispersal_params <- as.dispersal(list(dispersal_distance=list('larvae'=3,'juvenile'=0,'adult'=3),
#'                dispersal_kernel=list('larvae'=exp(-c(0:2)),'juvenile'=0,'adult'=exp(-c(0:2)*.2)),
#'                dispersal_proportion=list('larvae'=0.1,'juvenile'=0,'adult'=0.3)))  
#'
#' ## dispersal using cellular automata.                                 
#' dispersed_populations <- dispersal(dispersal_params, habitat, method='ca')
#' populations(habitat) <- dispersed_populations
#' 
#' ## dispersal using fast fourier transforms (diffuse).
#' dispersed_populations <- dispersal(dispersal_params, habitat, method='fft')
#' populations(habitat) <- dispersed_populations


as.dispersal <- function (params) {
  stopifnot(is.list(params))
  stopifnot(length(params)>=3)
  if(!exists(c('dispersal_distance','dispersal_kernel','dispersal_proportion'),params))stop('dispersal parameters must contain "dispersal_distance","dispersal_kernel" and "dispersal_proportion"')
  class(params) <- 'dispersal'
  return(params)
}

#' @rdname dispersal-class
#' @export
#' @examples
#' is.dispersal(disp)
is.dispersal <- function (x) inherits(x, 'dispersal')

#' @rdname dispersal-class
#' @param x an object to be tested as a dispersal transfun object
#' @param \dots further arguments passed to or from other methods.
#' @export
#' @examples
#' # print method
#' print(dp)
#'
print.dispersal <- function(x,...){
  
  vals <- which(names(x)%in%c("dispersal_distance","dispersal_kernel","dispersal_proportion"))
  x <- x[vals]
  
  disp_info <- list()
  for (i in base::seq_along(x)) {
  disp_info[[i]]  <- base::paste(
      base::names(x[[i]]),
      base::format(x[[i]], digits= 3,scientific = FALSE),
      sep = " = ",
      collapse = "\n "
    )
  };
  text <- sprintf('dispersal function with disperal distance of:\n %s \ncells;\n\ndispersal kernels of:\n %s;\n\nand a dispersal proportion of:\n %s\n for stages.',disp_info[[1]],disp_info[[2]],disp_info[[3]])
  cat(text)
}

#' @rdname dispersal-class
#' @name extent
#' @export

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

#' @rdname dispersal-class
#' @name bcb
#' @export

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

#' @rdname dispersal-class
#' @name setupFFT
#' @export

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

#' @rdname dispersal-class
#' @name dispersalFFT
#' @export
#' @examples 
#' # small functions that should be in base R
#'seq_range <- function (range, by = 1) seq(range[1], range[2], by = by)
#'ifft <- function (z) fft(z, inverse = TRUE)
#'
#'# coordinates of grid cells
#'# *don't make n too big, as the dense version will take forever!*
#'n <- 50 * c(2, 2)
#'y <- seq_len(n[1])
#'x <- seq_len(n[2])
#'
#'# dispersal function acting on distance matrix
#'# cut-off dispersal at the minimum dimension of the grid 
#'f <- function (d, cutoff = min(n)) {
#'  ifelse (d > cutoff, 0, exp(-d))
#'}
#'
#'# f <- function (d) exp(-d)
#'# initial population on grid (one stage)
#'pop <- matrix(rpois(length(x) * length(y),10),
#'              length(y), length(x))
#'              
#'# setup for the fft approach (run this once, before the simulation)
#'fs <- setupFFT(x = x, y = y, f = f, factor = 1)
#'
#'# apply dispersal to the population (need to run this separately for each stage)
#'pop_new <- dispersalFFT(popmat = pop, fs = fs)

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
  pop_fft <- fft(t(fs$pop_torus))
  bcb_fft <- fft(fs$bcb_vec)
  pop_new_torus_fft <- ifft(pop_fft * bcb_fft)
  
  # convert back to real domain, apply correction and transpose
  pop_torus_new <- t(Re(pop_new_torus_fft / length(fs$pop_torus)))
  
  # extract the section of the torus representing our 2D plane and return
  pop_new <- pop_torus_new[fs$yidx, fs$xidx]
  pop_new
}

#' @rdname dispersal-class
#' @export
seq_range <- function (range, by = 1) seq(range[1], range[2], by = by)

#' @rdname dispersal-class
#' @export
ifft <- function (z) fft(z, inverse = TRUE)

#' @rdname dispersal-class
#' @name dispersal
#' @export 

dispersal <- function(params,habitat,method,...){
                          stopifnot(is.list(params))
                          if(!any(method==c('ca','fft')))stop('method must be either "ca" (cellular automata) or\n "fft" (fast fourier transformation).')
                          stopifnot(is.habitat(habitat))  
                          dispersal_results <- switch(method,
                                                     ca = dispersal_core_ca(params,habitat),
                                                     fft = dispersal_core_fft(params,habitat))  
                          return(dispersal_results)
}

dispersal_core_ca <- function(params,habitat){

  #generate default parameters for dispersal parameters if they are missing from 'params'. 
  if(!exists('barrier_type',params))params$barrier_type <- 0
  if(!exists('dispersal_steps',params))params$dispersal_steps <- 1
  if(!exists('use_barriers',params))params$use_barriers <- FALSE
  
  #identify populations and workout which populations can disperse.
  which_stages_disperse <- which(params$dispersal_proportion>0)
  n_dispersing_stages <- length(which_stages_disperse)
  
  #extract habitat suitability, relevant populations and carrying capacity. 
  pops <- populations(habitat)
  disperse_pops <- pops[which_stages_disperse]
  cc <- carrying_capacity(habitat)
  hsm <- habitat_suitability(habitat)
  
  #if barriers is NULL create a barriers matrix all == 0.
  if(!exists('barriers_map',params)){
    bm <- raster::calc(hsm,function(x){x[!is.na(x)] <- 0; return(x)})
    params$barriers_map <- bm
  }
  
  ca_dispersal <- list()
  # could do this in parallel if wanted. 
  for (i in seq_len(n_dispersing_stages)){
                              ca_dispersal[[i]] <- dhmpr::rcpp_dispersal(raster::as.matrix(disperse_pops[[i]]), raster::as.matrix(cc), raster::as.matrix(hsm),raster::as.matrix(params$barriers_map), as.integer(params$barrier_type), params$use_barrier, as.integer(params$dispersal_steps), as.integer(params$dispersal_distance[which_stages_disperse][i]), as.numeric(unlist(params$dispersal_kernel[which_stages_disperse][i])), as.numeric(params$dispersal_proportion[which_stages_disperse][i]))[[1]] # we only want the dispersal population matricies.
  }
  
  ca_dispersal <- lapply(ca_dispersal,function(x){hsm[]<-x;return(hsm)})
  pops[which_stages_disperse] <- ca_dispersal
  pops <- lapply(pops, `attr<-`, "habitat", "populations")
  return(pops)
}

#### up to here <----
dispersal_core_fft <- function(params,habitat){
   
  #identify populations and workout which populations can disperse.
  which_stages_disperse <- which(params$dispersal_proportion>0)
  n_dispersing_stages <- length(which_stages_disperse)
  
  ## get the relevant 
  pops <- populations(habitat)
  disperse_pops <- pops[which_stages_disperse]
  n <- dim(pops[[1]])[1:2]
  y <- seq_len(n[1])
  x <- seq_len(n[2])
  
  ## set up disperal function
  f <- function (d, cutoff = min(n)) {
                   ifelse (d > cutoff, 0, exp(-d))
                  }
                  
  # f <- function (d) exp(-d)
  # setup for the fft approach (run this once, before the simulation)
  fs <- setupFFT(x = x, y = y, f = f)
                  
  #'# apply dispersal to the population (need to run this separately for each stage)
  fft_dispersal <- list()
  # could do this in parallel if wanted. 
  for (i in seq_len(n_dispersing_stages)){
    fft_dispersal[[i]] <- dispersalFFT(popmat = raster::as.matrix(pops[[i]]), fs = fs)
  }
  
  fft_dispersal <- lapply(fft_dispersal,function(x){pops[[1]][]<-x;return(pops[[1]])})
  pops[which_stages_disperse] <- fft_dispersal
  pops <- lapply(pops, `attr<-`, "habitat", "populations")
  return(pops)
}
