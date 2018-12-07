#' Run an simulation
#'
#' A simulation changes landscape objects based on selected dynamics over a
#' specified number of timesteps.
#'
#' @rdname simulation_results
#'
#' @param landscape a landscape object representing the initial habitat and
#'   population
#' @param population_dynamics a population_dynamics object describing how
#'   population changes over time
#' @param habitat_dynamics optional list of functions to modify the landscape at
#'   each timestep
#' @param timesteps number of timesteps used in one simulation or to display
#'   when plotting rasters
#' @param replicates number simulations to perform
#' @param verbose print messages and progress to console? (default is TRUE)
#' @param x an simulation_results object
#' @param object the simulation_results object to plot - can be 'population'
#'   (default), 'suitability' or 'carrying_capacity'
#' @param type the plot type - 'graph' (default) or 'raster'
#' @param stages life-stages to plot - must be specified for 'raster' plot types;
#'   default is NULL and all life-stages will be plotted
#' @param animate if plotting type 'raster' would you like to animate the
#'   rasters as a gif?
#' @param panels the number of columns and rows to use when plotting raster
#'   timeseries - default is 3 x 3 (e.g. c(3,3) )
#' @param replicate which replicate to extract a spatial object from a simulation result
#' @param timestep which timestep to extract a spatial object from a simulation result
#' @param landscape_object which landscape object to extract a spatial object from a
#'   simulation result - can be specified by name (e.g. "suitability") or index number
#' @param stage which life-stage to extract a spatial object from a simulation result
#'   - only used for 'population' components of the landscape object
#' @param misc which misc object to extract a spatial object from a simulation result
#'   - only used for additional spatial objects added to a landscape (e.g. disturbance layers)
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{simulation_results}
#'
#' @export
#' 
#' @importFrom future plan multiprocess future values
#' @importFrom raster animate
#' @importFrom viridisLite viridis
#'
#' @examples
#'
#' library(steps)
#' library(raster)
#' library(future)
#' plan(multiprocess)
#'
#' landscape <- landscape(population = egk_pop,
#'                suitability = egk_hab, 
#'                carrying_capacity = egk_k)
#'
#' pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat))
#'
#' results <- simulation(landscape, pop_dynamics, timesteps = 10, replicates = 5)

simulation <- function(landscape, population_dynamics, habitat_dynamics = list(), timesteps = 3, replicates = 1, verbose = TRUE){

  in_parallel <- !inherits(future::plan(), "sequential")
  lapply_fun <- ifelse(in_parallel,
                       future.apply::future_lapply,
                       base::lapply)
  
  simulation_results <- lapply_fun(seq_len(replicates),
                                      FUN = simulate,
                                      landscape = landscape,
                                      population_dynamics = population_dynamics,
                                      habitat_dynamics = habitat_dynamics,
                                      timesteps = timesteps,
                                      verbose = verbose)

  as.simulation_results(simulation_results)
}


#' @export
#' @noRd
`[.simulation_results` <- function(x, ..., drop=TRUE) {
  structure(NextMethod(), class=class(x))
}


#' @rdname simulation_results
#'
#' @export
#'
#' @examples
#'
#' # Test if object is of the type 'simulation_results'
#'
#' is.simulation_results(results)

is.simulation_results <- function (x) {
  inherits(x, 'simulation_results')
}

#' @rdname simulation_results
#'
#' @export
#'
#' @examples
#'
#' print(results)

print.simulation_results <- function (x, ...) {
  cat("This is an simulation results object, for", length(x), "replicates")
  
}


#' @rdname simulation_results
#'
#' @export
#'
#' @examples
#'
#' plot(results)

plot.simulation_results <- function (x,
                                     object = "population",
                                     type = "graph",
                                     stages = NULL,
                                     animate = FALSE,
                                     timesteps = c(1:9),
                                     panels = c(3,3),
                                     ...){
  
  total_stages <- raster::nlayers(x[[1]][[1]]$population)
  stage_names <- names(x[[1]][[1]]$population)

  graph.pal <- c(
    "#b4d7b3",
    "#eb7d75",
    "#80b1d3",
    "#bebada",
    "#f0ab7e",
    "#969696"
  )
  
  if (length(x) == 1) {
   
    if (object == "population") {
      
      pop <- get_pop_replicate(x[[1]])
      
      if (type == "graph") {
        
        if (is.null(stages)) {
          
          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1, total_stages))
          
          for (i in seq_len(total_stages)) {
            
            graphics::plot(pop[ , i],
                           type='l',
                           ylab=paste("Total Population: ",stage_names[i]),
                           xlab="Timesteps",
                           lwd=3,
                           col=graph.pal[i],
                           ylim=range(pretty(pop)))

          }
          
        }
        
        if (!is.null(stages) && stages == 0) {
          
          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
          
          graphics::plot(rowSums(pop),
                         type='l',
                         ylab="Total Population (all stages)",
                         xlab="Timesteps",
                         lwd=3,
                         col="black",
                         ylim=range(pretty(rowSums(pop))))

        }
        
        if (!is.null(stages) && length(stages) > 0) {
          
          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,length(stages)))
          
          for (i in stages) {
          
            graphics::plot(pop[, i],
                           type='l',
                           ylab=paste("Total Population: ",stage_names[stages]),
                           xlab="Timesteps",
                           lwd=3,
                           col=graph.pal[i],
                           ylim=range(pretty(pop)))
          }

        }
        
      }
      
      if (type == "raster") {
        
        if (is.null(stages)) {
          stop("Please provide a life-stage when plotting population rasters or specify zero (0) for a sum of all life-stages")
        }
        
        if (length(timesteps) > 20) {
          cat ("Note, you have specified to plot rasters for more than 20 timesteps;",
               "\nthis could take a while depending on the resolution of your landscape.",
               "\nPress [enter] to continue, or [esc] to cancel and adjust settings.")
          line <- readline()
        }
        
        if(stages == 0) {
          
          rasters_sum <- raster::stack(lapply(x[[1]], function (landscape) sum(landscape$population)))
          
          # Find maximum and minimum population value in raster cells for all timesteps for life-stage
          scale_max <- ceiling(max(raster::cellStats(rasters_sum, max)))
          scale_min <- floor(min(raster::cellStats(rasters_sum, min)))
          
          # Produce scale of values
          breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)
          
          if (animate == TRUE) {
            graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
            
            raster::animate(rasters_sum[[timesteps]],
                            col=viridisLite::viridis(length(breaks)-1),
                            n = 1)
          } else {
            
            graphics::par(mar=c(2, 0, 0, 0), mfrow=c(1,1))
            print(rasterVis::levelplot(rasters_sum[[timesteps]],
                                       scales = list(draw = FALSE),
                                       margin = list(draw = FALSE),
                                       at = breaks,
                                       col.regions = viridisLite::viridis(length(breaks)-1),
                                       colorkey = list(space = "bottom",
                                                       title = "individuals",
                                                       width = 0.4),
                                       main = "population",
                                       layout = panels))
          }
          
        } else {
          
          rasters <- raster::stack(lapply(x[[1]], function (landscape) landscape$population[[stages]]))
          
          # Find maximum and minimum population value in raster cells for all timesteps for life-stage
          scale_max <- ceiling(max(raster::cellStats(rasters, max)))
          scale_min <- floor(min(raster::cellStats(rasters, min)))
          
          # Produce scale of values
          breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)
          
          if (animate == TRUE) {
            graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
            
            raster::animate(rasters[[timesteps]],
                            col=viridisLite::viridis(length(breaks)-1),
                            n = 1)
          } else {

            print(rasterVis::levelplot(rasters[[timesteps]],
                                       scales = list(draw = FALSE),
                                       margin = list(draw = FALSE),
                                       at = breaks,
                                       col.regions = viridisLite::viridis(length(breaks)-1),
                                       colorkey = list(space = "bottom",
                                                       title = "individuals",
                                                       width = 0.4),
                                       main = "population",
                                       layout = panels))
            
          }
        }
      }
    }
    
    if (object == "suitability") {
      
      rasters <- raster::stack(lapply(x[[1]], function (landscape) landscape$suitability))
      
      # Find maximum and minimum population value in raster cells for all timesteps for life-stage
      scale_max <- ceiling(max(raster::cellStats(rasters, max)))
      scale_min <- floor(min(raster::cellStats(rasters, min)))
      
      # Produce scale of values
      breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)
      
      if (animate == TRUE) {
        graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
        
        raster::animate(rasters[[timesteps]],
                        col=viridisLite::viridis(length(breaks)-1),
                        n = 1)
      
      } else {  

        print(rasterVis::levelplot(rasters[[timesteps]],
                                   scales = list(draw = FALSE),
                                   margin = list(draw = FALSE),
                                   at = breaks,
                                   col.regions = viridisLite::viridis(length(breaks)-1),
                                   colorkey = list(space = "bottom",
                                                   title = "index",
                                                   width = 0.4),
                                   main = "habitat",
                                   layout = panels))
      }
        
    }
    
    if (object == "carrying_capacity") {
      
      rasters <- raster::stack(lapply(x[[1]], function (landscape) landscape$carrying_capacity))
      
      # Find maximum and minimum population value in raster cells for all timesteps for life-stage
      scale_max <- ceiling(max(raster::cellStats(rasters, max)))
      scale_min <- floor(min(raster::cellStats(rasters, min)))
      
      # Produce scale of values10
      breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)
      
      if (animate == TRUE) {
        graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
        
        raster::animate(rasters[[timesteps]],
                        col=viridisLite::viridis(length(breaks)-1),
                        n = 1)
        
      } else {  

        print(rasterVis::levelplot(rasters[[timesteps]],
                                   scales = list(draw = FALSE),
                                   margin = list(draw = FALSE),
                                   at = breaks,
                                   col.regions = viridisLite::viridis(length(breaks)-1),
                                   colorkey = list(space = "bottom",
                                                   title = "individuals",
                                                   width = 0.4),
                                   main = "k",
                                   layout = panels))
      }
        
    }
     
  }
  
  if (length(x) > 1) {
    
    pop <- get_pop_simulation(x)
    pop.mn <- round(apply(pop, c(1,2), mean), 0)
    
    if (type == "raster" | object == "suitability" | object == "carrying_capacity") {
      stop("Raster plotting is only available for single simulations")
    }
    
    if (is.null(stages)) {
      
      graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,total_stages))
      
      for (i in seq_len(total_stages)) {
        
        graphics::plot(pop.mn[, i],
                       type = 'l',
                       ylab = paste("Total Population: ",stage_names[i]),
                       xlab = "Timesteps",
                       lwd = 3,
                       col = graph.pal[i],
                       ylim=range(pretty(pop)))

        for (j in seq_along(x)) {
          graphics::lines(pop[ , i, j],
                          col = 'gray')
        }
        
      }
      
    }
    
    if(!is.null(stages) && stages == 0) {
      
      graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
      
      # draw the 95% CI polygon (if available) and median line
      quants <- t(apply(apply(pop, 3, rowSums),1, stats::quantile, c(0.025, 0.5, 0.975)))
      
      xaxs <- seq_len(nrow(pop[,,1]))
      
      graphics::plot(quants[,2],#rowSums(pop[ , , 1]),
                     type = 'l',
                     ylab = "Total Population (all stages)",
                     xlab = "Timesteps",
                     lwd = 3,
                     col = 'black')
      
      # for (j in seq_along(x)[-1]) {
      #   graphics::lines(rowSums(pop[ , , j]),
      #                   col = 'gray')
      # }
      
      graphics::polygon(x = c(xaxs, rev(xaxs)),
              y = c(quants[, 1], rev(quants[, 3])),
              col = grDevices::grey(0.9),
              border = NA)
      
      graphics::lines(quants[, 2] ~ xaxs,
            lwd = 2,
            col = grDevices::grey(0.4))
      
      if (!is.null(x[[1]][[1]]$carrying_capacity)) {
        graphics::abline(h=raster::cellStats(x[[1]][[1]]$carrying_capacity,sum),
                         lwd=1,
                         lty=2)
      }

    }
    
    if (!is.null(stages) && stages > 0) {
      
      graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
      
      graphics::plot(pop.mn[ , stages],
                     type = 'l',
                     ylab = paste("Total Population: ",stage_names[stages]),
                     xlab = "Timesteps",
                     lwd = 3,
                     col = graph.pal[stages])
      
      for (j in seq_along(x)) {
        graphics::lines(pop[ , stages, j],
                        col = 'gray')
      }

    }
  
  }
  
}

#' @rdname simulation_results
#'
#' @export
#'
#' @examples
#'
#' # Extract spatial components from a 'simulation_results' object
#'
#' extract_results(results)

extract_results <- function (x,
                            replicate = 1,
                            timestep = 1,
                            landscape_object = "population",
                            stage = 1,
                            misc = 1) {
  if (landscape_object == 1 | landscape_object == "population") {
    x[[replicate]][[timestep]][[landscape_object]][[stage]]
  }
  if (landscape_object > 3 |
      !landscape_object == "population" |
      !landscape_object == "suitability" |
      !landscape_object == "carrying_capacity") {
    x[[replicate]][[timestep]][[landscape_object]][[misc]]
  }

}

##########################
### internal functions ###
##########################

as.simulation_results <- function (simulation_results) {
  as_class(simulation_results, "simulation_results", "list")
}

simulate <- function (i, landscape, population_dynamics, habitat_dynamics, timesteps, verbose) {
  timesteps <- seq_len(timesteps)
  if (verbose == TRUE) pb <- utils::txtProgressBar(min = 0, max = max(timesteps), style = 3)
  
  output_landscapes <- list()
  
  for (timestep in timesteps) {
    
    landscape <- population_dynamics(landscape, timestep)
    
    for (dynamic_function in habitat_dynamics) {
      landscape <- dynamic_function(landscape, timestep)
    }
    
    output_landscapes[[timestep]] <- landscape
    if (verbose == TRUE) utils::setTxtProgressBar(pb, timestep)
  }
  
  if (verbose == TRUE) close(pb)
  
  as_class(output_landscapes, "replicate", "list")
}

# extract populations from a simulation
get_pop_replicate <- function(x, ...) {
  total_stages <- raster::nlayers(x[[1]]$population)
  idx <- which(!is.na(raster::getValues(x[[1]]$population[[1]])))
  pops <- lapply(x, function(x) raster::extract(x$population, idx))
  pop_sums <- lapply(pops, function(x) colSums(x))
  pop_matrix <- matrix(unlist(pop_sums), ncol = total_stages, byrow = TRUE)
  return(pop_matrix)
}

get_pop_simulation <- function(x, ...) {
  total_stages <- raster::nlayers(x[[1]][[1]]$population)
  timesteps <- length(x[[1]])
  sims <- length(x)
  
  pop_array <- array(dim=c(timesteps,total_stages,sims))
  
  for(i in seq_len(sims)) {
    pop_array[, , i] <- get_pop_replicate(x[[i]])
  }
  return(pop_array)
}
