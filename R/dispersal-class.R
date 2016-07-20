#' @title dispersal class for metapoplation
#' @rdname dispersal 
#' @name dispersal
#' @param params list which containts alpha the exponential decay rate of patch connectivity 
#' (dispersion parameter), beta double parameter that represents the shape of the dispersal kernel.
#' and disp_fun a characture 'H' uses hanski(1994), if 'S' uses shaw(1995).
#' @export
#' @examples 
#' params <- list(alpha=1,beta=1,disp_fun="H")
#' dispersal(params)
#' 
#' dispersal(NULL)

dispersal <- function(params){
  switch(class(params[1]),
         NULL = dispersalDefault(),
         list = dispersal_core(params))
}

#' @rdname dispersal
#' @name d
#' @export
#' @examples 
#' d(NULL)
#' 
d <- dispersal

dispersal_core <- function(params)
  {
  
  dist <- function(habitat) {
           dist <- distance(habitat)
           dist
  }
  params$dist <- dist
  return(params)
}


dispersalDefault <- function () {
  params_list <- list(alpha = 1,
                      beta = 1,
                      disp_fun = 'H')
  d <- dispersal(params_list)
  return (d)
}

#' @rdname dispersal
#' @param x an object to be tested as a dispersal transfun object
#' @export
#' @examples
#' is.dispersal(disp)
is.dispersal <- function (x) inherits(x, 'dispersal')

#' @rdname dispersal
#' @param x an object to be tested as a dispersal transfun object
#' @export
as.dispersal <- function (x) {
  if (!is.transition(x)) {
      class(x) <- c('transition', class(x))
    }
    if (!is.dispersal(x)) {
      class(x) <- c('dispersal', class(x))
    }
    return (x)
}
