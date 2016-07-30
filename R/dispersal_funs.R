#' @title dispersal functions
#' @rdname disp_funs
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

#' @rdname disp_funs
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

#' @rdname disp_funs
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

#' @rdname disp_funs
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