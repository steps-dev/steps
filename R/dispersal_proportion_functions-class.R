#' Create a proportion dispersing function
#'
#' @description A proportion dispersing function generates the proportions of species
#' that disperse from cells based on landscape features.
#'  
#' The default proportion dispersing function returns full dispersal (1.0) for all life stages.
#' Additional proportion dispersing functions are also provided in the software for the
#' user to select (see' \link[steps]{carrying_capacity_dispersal}), however, a user may also provide a
#' custom written proportion dispersing function.
#' 
#' @rdname dispersal_proportion_function
#'
#' @param proportions proportions of individuals in each life stage that disperse - default is 1
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
    
    dispersal_proportion <- rep(1, n_stages)
    
    if (length(proportions) > n_stages | length(proportions) < n_stages)  {
      if (timestep == 1) cat("    ",
                             n_stages,
                             "life stages exist but",
                             length(proportions),
                             "dispersal proportion(s) of",
                             paste(proportions, collapse = ", "),
                             "were specified. Please check your values.")
      dispersal_proportion <- rep_len(proportions, n_stages)
    } else {
      dispersal_proportion <- proportions
    }
    
    dispersal_proportion
    
  }
  
  as.dispersal_proportion_function(disp_prop_fun)
}


#' @rdname dispersal_proportion_function
#' 
#' @export
#'
#' @examples
#' 
#' test_dispersal_function <- carrying_capacity_dispersal()

carrying_capacity_dispersal <- function () {
  
  disp_prop_fun <- function (landscape, timestep) {
    
    # get total life-stages
    n_stages <- raster::nlayers(landscape$population)
    
    # check for specified stages that contribute to density dependence
    if (!exists("density_dependence_stages")) {
      density_dependence_stages <- seq_len(n_stages)
    }
    
    # get non-NA cells
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))

    pop <- raster::getValues(landscape$population)
    
    cc <- get_carrying_capacity(landscape, timestep)
    cc <- raster::getValues(cc)
    
    dispersal_proportion <- rep(0, n_stages)
    
    for(i in density_dependence_stages) {
      proportion <- sum(pop[cell_idx, i]) / sum(cc[cell_idx])
      dispersal_proportion[i] <- tanh(proportion)
    }
    
    dispersal_proportion
    
  }
  
  as.dispersal_proportion_function(disp_prop_fun)
}


##########################
### internal functions ###
##########################

as.dispersal_proportion_function <- function (dispersal_proportion_function) {
  as_class(dispersal_proportion_function, "dispersal_proportion_function", "function")
}
