#' Create a proportion dispersing function
#'
#' @description A proportion dispersing function generates the proportions of species
#' that disperse from cells based on landscape features.
#'  
#' The default proportion dispersing function returns full dispersal (1.0) for all life stages.
#' Additional proportion dispersing functions are also provided in the software for the
#' user to select (see' \link[steps]{exponential_dispersal_kernel}), however, a user may also provide a
#' custom written dispersal kernel.
#' 
#' @rdname dispersal_proportion_function
#'
#' @param proportions which controls the rate at which the population disperses with distance
#' @param normalize (exponential dispersal parameter) should the normalising constant be used - default is FALSE.
#' 
#' @return An object of class \code{dispersal_proportion_function}
#' 
#' @export
#'
#' @examples
#' 
#' test_dispersal_function <- all_dispersing()

all_dispersing <- function (proportions = 1) {
  
  disp_prop_fun <- function (landscape, timestep) {
   
    n_stages <- raster::nlayers(landscape$population)
    
    if (length(proportions) > n_stages) {
      if (timestep == 1) cat("    ",
                             n_stages,
                             "life stages exist but",
                             length(dispersal_proportion),
                             "dispersal proportion(s) of",
                             dispersal_proportion,
                             "were specified. Please check your values")
      dispersal_proportion <- rep(proportions, n_stages)
    } else {
      dispersal_proportion <- proportions
    }
    
    dispersal_proportion
    
  }
  
  as.dispersal_proportion_function(disp_prop_fun)
}


get_proportion_dispersing <- function(proportion_dispersing, landscape, timestep) {
  
  pd <- proportion_dispersing
  
  if (inherits(pd, "numeric")) {
    
    res <- cc
    
  } else if (is.function(cc)) {
    
    # if it's a function, run it on landscape
    res <- pd(landscape, timestep)
    
  } else {
    
    # otherwise, we don't support it
    stop ("invalid proportion dispersing argument",
          call. = FALSE)
    
  }
  
  res
  
}


##########################
### internal functions ###
##########################

as.dispersal_proportion_function <- function (dispersal_proportion_function) {
  as_class(dispersal_proportion_function, "dispersal_proportion_function", "function")
}
