#' Run a simulation
#'
#' A simulation changes landscape objects based on selected dynamics over a
#' specified number of timesteps.
#'
#' @rdname simulation
#'
#' @param landscape a \link[steps]{landscape} object representing the initial habitat and
#'   population
#' @param population_dynamics a \link[steps]{population_dynamics} object describing how
#'   population changes over time
#' @param habitat_dynamics optional list of functions to modify the landscape at
#'   each timestep - see \link[steps]{habitat_dynamics_functions}
#' @param demo_stochasticity how should population rounding occur, if at all -
#'   "full" uses a multinomial draw to return rounded cell populations (default)
#'   whilst "none" returns non-integer cell populations (no rounding). Note, this
#'   parameter specification is used consistently throughout all functions in a
#'   simulation.
#' @param timesteps number of timesteps used in one simulation
#' @param replicates number of simulations to perform
#' @param verbose print messages and progress to console? (default is TRUE)
#' @param future.globals a list of custom functions, and objects called by the functions,
#'   that a user has created in the global environment for use in a simulation. Note this
#'   is only required when running simulations in parallel (e.g. plan(multisession)). 
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
#' \dontrun{
#' ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)
#' 
#' pd <- population_dynamics(change = growth(egk_mat),
#'                           dispersal = kernel_dispersal(max_distance = 2000,
#'                                         dispersal_kernel = exponential_dispersal_kernel(
#'                                           distance_decay = 1000)),
#'                           density_dependence = ceiling_density())
#' 
#' # Run a simulation with full demographic stochasticity and without any habitat
#' # dynamics for tewnty timesteps.
#' sim <- simulation(landscape = ls,
#'                   population_dynamics = pd,
#'                   habitat_dynamics = NULL,
#'                   timesteps = 20)
#' }

simulation <- function(landscape,
                       population_dynamics,
                       habitat_dynamics = list(),
                       demo_stochasticity = c("full", "none"),
                       timesteps = 3,
                       replicates = 1,
                       verbose = TRUE,
                       future.globals = list()){
  
  # gather globals in future.apply call
  future_lapply_wrapper <- function (...) {
    future.apply::future_lapply(..., future.globals = future.globals)
  }
  
  # clear out the stash every time we begin a simulation
  flush_stash()
  
  steps_stash$demo_stochasticity <- match.arg(demo_stochasticity)
  
  in_parallel <- !inherits(future::plan(), "sequential")
  is_multisession <- inherits(future::plan(), "multisession")
  lapply_fun <- ifelse(in_parallel,
                       future_lapply_wrapper,
                       base::lapply)
  
  landscape_names <- names(landscape)
  
  # loop through names of objects and check for layer counts
  for (name in landscape_names) {
    
    if (name != "population" && !is.function(landscape[[name]]) &&
        !is.null(landscape[[name]]) && raster::nlayers(landscape[[name]]) > 1) {
      
      if (raster::nlayers(landscape[[name]]) < timesteps) {
        stop("A spatial object exists in the landscape that has less layers than specified ",
             "number timesteps. All spatial objects must have either one layer or a number of ",
             "layers equal to the intended number of timesteps in a simulation. Please check the ",
             "landscape object and re-run the simulation.")
      }
      
    }
    
  }
  
  simulation_results <- tryCatch(lapply_fun(seq_len(replicates),
                                            FUN = simulate,
                                            landscape = landscape,
                                            population_dynamics = population_dynamics,
                                            habitat_dynamics = habitat_dynamics,
                                            timesteps = timesteps,
                                            verbose = verbose,
                                            stash = steps_stash,
                                            is_multisession = is_multisession),
                                 error = global_object_error)
  
  as.simulation_results(simulation_results)
}

#' @export
#' @noRd
`[.simulation_results` <- function(x, ..., drop=TRUE) {
  structure(NextMethod(), class=class(x))
}

is.simulation_results <- function (x) {
  inherits(x, 'simulation_results')
}

# #' @rdname simulation
# #'
# #' @export
# #'
# #' @examples
# #'
# #' print(results)
# 
# print.simulation_results <- function (x, ...) {
# 
#   cat("This is an simulation results object, for", length(x), "replicates")
# 
# }


#' Plot the results of a simulation
#' 
#' Methods to visually inspect the results of a simulation. Both linear graphs
#' and spatial-explicit grids can be generated to illustrate population changes through
#' time and space.
#' 
#' @param x a simulation_results object
#' @param object the \code{simulation_results} object to plot - can be 'population'
#'   (default), 'suitability' or 'carrying_capacity'
#' @param type the plot type - 'graph' (default) or 'raster'
#' @param stages life-stages to plot - by default all life-stages will be considered.
#'   Set to zero for totals (i.e. sum of all life-stages). For raster plotting,
#'   the life-stages that are specified will be summed, unless a single life-stage
#'   is specified.
#' @param animate if plotting type 'raster' would you like to animate the
#'   rasters?
#' @param timesteps timesteps to display when plotting rasters
#' @param panels the number of columns and rows to use when plotting raster
#'   timeseries - default is 3 x 3 (e.g. c(3, 3) )
#' @param emp should the expected minimum population of the simulation be
#'   plotted?
#' @param ... further arguments passed to/from other methods
#'   
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)
#' 
#' pd <- population_dynamics(change = growth(egk_mat),
#'                           dispersal = kernel_dispersal(max_distance = 2000,
#'                                         dispersal_kernel = exponential_dispersal_kernel(
#'                                           distance_decay = 1000)),
#'                           density_dependence = ceiling_density())
#' 
#' sim <- simulation(landscape = ls,
#'                   population_dynamics = pd,
#'                   habitat_dynamics = NULL,
#'                   timesteps = 20)
#' 
#' # Plot the population trajectories by life-stage
#' plot(sim)
#' 
#' # Plot the spatial distributions of total populations for first nine timesteps
#' plot(sim, type = "raster", stages = 0, timesteps = 1:9) 
#' }

plot.simulation_results <- function (x,
                                     object = "population",
                                     type = "graph",
                                     stages = NULL,
                                     animate = FALSE,
                                     timesteps = c(1:3),
                                     panels = c(3,3),
                                     emp = FALSE,
                                     ...){
  
  # don't have a persistent effect on the graphics device
  # op <- graphics::par(no.readonly = TRUE)
  # on.exit(graphics::par(op))
  
  total_stages <- raster::nlayers(x[[1]][[1]]$population)
  stage_names <- names(x[[1]][[1]]$population)
  
  graph.pal <- c(
    "#6da36b",
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
                           type = 'l',
                           ylab = paste("Total Population: ", stage_names[i]),
                           xlab = "Timesteps",
                           #lwd = 3,
                           col = graph.pal[i],
                           ylim = range(pretty(pop)),
                           xaxt = 'n',
                           xlim = c(1, length(pop[ , i])),
                           ...)
            
            graphics::axis(side = 1, at = pretty_int(c(1, length(pop[ , i]))))
            
          }
          
        }
        
        if (!is.null(stages) && stages == 0) {
          
          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1))
          
          graphics::plot(rowSums(pop),
                         type='l',
                         ylab="Total Population (all stages)",
                         xlab="Timesteps",
                         #lwd=3,
                         col="black",
                         xaxt = 'n',
                         xlim = c(1, length(rowSums(pop))),
                         ...)
          
          graphics::axis(side = 1, at = pretty_int(c(1, length(rowSums(pop)))))
          
        }
        
        if (!is.null(stages) && length(stages) > 0 && stages != 0) {
          
          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1, length(stages)))
          
          for (i in stages) {
            
            graphics::plot(pop[ , i],
                           type='l',
                           ylab=paste("Total Population: ",stage_names[stages]),
                           xlab="Timesteps",
                           #lwd=3,
                           col=graph.pal[i],
                           ylim=range(pretty(pop)),
                           xaxt = 'n',
                           xlim = c(1, length(pop[ , i])),
                           ...)
            
            graphics::axis(side = 1, at = pretty_int(c(1, length(pop[ , i]))))
            
          }
          
        }
        
      }
      
      if (type == "raster") {
        
        if(!is.null(stages) && stages == 0) {
          
          rasters_sum <- raster::stack(lapply(x[[1]], function (landscape) sum(landscape$population)))
          #rasters_sum[rasters_sum == 0] <- NA
          
          names(rasters_sum) <- paste0("Timestep_", 1:raster::nlayers(rasters_sum))
          
          # Find maximum and minimum population value in raster cells for all timesteps for life-stage
          scale_max <- ceiling(max(raster::cellStats(rasters_sum[[timesteps]], max)))
          scale_min <- floor(min(raster::cellStats(rasters_sum[[timesteps]], min)))
          #scale_max <- ceiling(max(stats::na.omit(raster::cellStats(rasters_sum, max))))
          #scale_min <- floor(min(stats::na.omit(raster::cellStats(rasters_sum, min))))
          
          # Produce scale of values
          breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)
          
          ifelse(any(rasters_sum[[timesteps]][] == 0, na.rm = TRUE),
                 colour_range <- c("#bfbfbfff", viridisLite::viridis(length(breaks)-1)),
                 colour_range <- viridisLite::viridis(length(breaks)-1))
          
          if (animate == TRUE) {
            graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
            
            raster::animate(rasters_sum[[timesteps]],
                            col = colour_range,
                            n = 1,
                            legend.args = list(text = 'individuals'),
                            box = FALSE,
                            axes = FALSE)
          } else {
            
            graphics::par(mar=c(2, 0, 0, 0), mfrow=c(1,1))
            print(rasterVis::levelplot(rasters_sum[[timesteps]],
                                       scales = list(draw = FALSE),
                                       margin = list(draw = FALSE),
                                       at = breaks,
                                       col.regions = colour_range,
                                       colorkey = list(space = "bottom",
                                                       width = 0.4),
                                       #main = "population",
                                       par.settings=list(layout.heights = list(xlab.key.padding = 1),
                                                         strip.background = list(col = "#ffffff")),
                                       xlab = "individuals",
                                       layout = panels))
          }
          
        } else {
          
          if (is.null(stages)) {
            stages <- seq_len(total_stages)
          }
          
          if (length(stages) == 1) {
            rasters <- raster::stack(lapply(x[[1]], function (landscape) landscape$population[[stages]]))
          } else {
            rasters <- raster::stack(lapply(x[[1]], function (landscape) sum(landscape$population[[stages]])))            
          }
          
          names(rasters) <- paste0("Timestep_", 1:raster::nlayers(rasters))
          
          # Find maximum and minimum population value in raster cells for all timesteps for life-stage
          scale_max <- ceiling(max(raster::cellStats(rasters[[timesteps]], max)))
          scale_min <- floor(min(raster::cellStats(rasters[[timesteps]], min)))
          #scale_max <- ceiling(max(stats::na.omit(raster::cellStats(rasters, max))))
          #scale_min <- floor(min(stats::na.omit(raster::cellStats(rasters, min))))
          
          # Produce scale of values
          breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)
          
          ifelse(any(rasters[[timesteps]][] == 0, na.rm = TRUE),
                 colour_range <- c("#bfbfbfff", viridisLite::viridis(length(breaks)-1)),
                 colour_range <- viridisLite::viridis(length(breaks)-1))
          
          if (animate == TRUE) {
            graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
            
            raster::animate(rasters[[timesteps]],
                            col = colour_range,
                            n = 1,
                            legend.args = list(text = 'individuals'),
                            box = FALSE,
                            axes = FALSE)
          } else {
            
            print(rasterVis::levelplot(rasters[[timesteps]],
                                       scales = list(draw = FALSE),
                                       margin = list(draw = FALSE),
                                       at = breaks,
                                       col.regions = colour_range,
                                       colorkey = list(space = "bottom",
                                                       width = 0.4),
                                       #main = "population",
                                       par.settings=list(layout.heights = list(xlab.key.padding = 1),
                                                         strip.background = list(col = "#ffffff")),
                                       xlab = "individuals",
                                       layout = panels))
            
          }
        }
      }
    }
    
    if (object == "suitability") {
      
      rasters <- raster::stack(lapply(x[[1]], function (x) x$suitability))
      
      names(rasters) <- paste0("Timestep_", 1:raster::nlayers(rasters))
      
      # Find maximum and minimum population value in raster cells for all timesteps for life-stage
      scale_max <- ceiling(max(stats::na.omit(raster::cellStats(rasters[[timesteps]], max))))
      scale_min <- floor(min(stats::na.omit(raster::cellStats(rasters[[timesteps]], min))))
      
      # Produce scale of values
      breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)
      
      ifelse(any(rasters[[timesteps]][] == 0, na.rm = TRUE),
             colour_range <- c("#bfbfbfff", viridisLite::viridis(length(breaks)-1)),
             colour_range <- viridisLite::viridis(length(breaks)-1))
      
      if (animate == TRUE) {
        graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
        
        raster::animate(rasters[[timesteps]],
                        col = colour_range,
                        n = 1,
                        legend.args = list(text = 'index'),
                        box = FALSE,
                        axes = FALSE)
        
      } else {  
        
        print(rasterVis::levelplot(rasters[[timesteps]],
                                   scales = list(draw = FALSE),
                                   margin = list(draw = FALSE),
                                   at = breaks,
                                   col.regions = colour_range,
                                   colorkey = list(space = "bottom",
                                                   width = 0.4),
                                   main = "habitat",
                                   par.settings=list(layout.heights = list(xlab.key.padding = 1),
                                                     strip.background = list(col = "#ffffff")),
                                   xlab = "suitability index",
                                   layout = panels))
      }
      
    }
    
    if (object == "carrying_capacity") {
      
      if (type == "graph") {
        
        idx <- which(!is.na(raster::getValues(x[[1]][[1]][["carrying_capacity"]])))
        k <- lapply(x[[1]], function(x) sum(raster::extract(x$carrying_capacity, idx)))
        
        graphics::par(mar=c(5.1, 4.1, 4.1, 2.1))
        
        graphics::plot(unlist(k),
                       type='l',
                       ylab="Total k (all stages)",
                       xlab="Timesteps",
                       #lwd=3,
                       col="black",
                       ylim=range(pretty(unlist(k))),
                       xaxt = 'n',
                       xlim = c(1, length(unlist(k))),
                       ...)
        
        graphics::axis(side = 1, at = pretty_int(c(1, length(unlist(k)))))
        
      }
      
      if (type == "raster") {
        
        rasters <- raster::stack(lapply(x[[1]], function (x) x$carrying_capacity))
        
        names(rasters) <- paste0("Timestep_", 1:raster::nlayers(rasters))
        
        # Find maximum and minimum population value in raster cells for all timesteps for life-stage
        scale_max <- ceiling(max(stats::na.omit(raster::cellStats(rasters[[timesteps]], max))))
        scale_min <- floor(min(stats::na.omit(raster::cellStats(rasters[[timesteps]], min))))
        
        # Produce scale of values10
        breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)
        
        ifelse(any(rasters[[timesteps]][] == 0, na.rm = TRUE),
               colour_range <- c("#bfbfbfff", viridisLite::viridis(length(breaks)-1)),
               colour_range <- viridisLite::viridis(length(breaks)-1))
        
        if (animate == TRUE) {
          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
          
          raster::animate(rasters[[timesteps]],
                          col = colour_range,
                          n = 1,
                          legend.args = list(text = 'individuals'),
                          box = FALSE,
                          axes = FALSE)
          
        } else {  
          
          print(rasterVis::levelplot(rasters[[timesteps]],
                                     scales = list(draw = FALSE),
                                     margin = list(draw = FALSE),
                                     at = breaks,
                                     col.regions = colour_range,
                                     colorkey = list(space = "bottom",
                                                     width = 0.4),
                                     main = "k",
                                     par.settings=list(layout.heights = list(xlab.key.padding = 1),
                                                       strip.background = list(col = "#ffffff")),
                                     xlab = "individuals",
                                     layout = panels))
        }
        
      }
    } 
  }
  
  if (length(x) > 1) {
    
    pop <- get_pop_simulation(x)
    pop.mn <- round(apply(pop, c(1,2), mean), 0)
    
    if (type == "raster") {
      stop("Raster plotting is only available for single replicates of simulations")
    }
    
    if (is.null(stages)) {
      
      graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1, total_stages))
      
      for (i in seq_len(total_stages)) {
        
        graphics::plot(pop.mn[, i],
                       type = 'l',
                       ylab = paste("Total Population: ", stage_names[i]),
                       xlab = "Timesteps",
                       #lwd = 3,
                       col = graph.pal[i],
                       ylim=range(pretty(pop)),
                       xaxt = 'n',
                       xlim = c(1, length(pop.mn[, i])),
                       ...)
        
        graphics::axis(side = 1, at = pretty_int(c(1, length(pop.mn[, i]))))
        
        for (j in seq_along(x)) {
          graphics::lines(pop[ , i, j],
                          col = 'gray',
                          lwd = 0.5)
        }
        
        graphics::lines(pop.mn[, i],
                        #lwd = 3,
                        col = graph.pal[i])
        
      }
      
    }
    
    if(!is.null(stages) && stages == 0) {
      
      graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
      
      # draw the 95% CI polygon (if available) and median line
      quants <- t(apply(apply(pop, 3, rowSums),1, stats::quantile, c(0.025, 0.5, 0.975)))
      
      xaxs <- seq_len(nrow(pop[ , , 1]))
      
      graphics::plot(quants[, 2], #rowSums(pop[ , , 1]),
                     type = 'l',
                     ylab = "Total Population (all stages)",
                     xlab = "Timesteps",
                     #lwd = 3,
                     col = 'black',
                     xaxt = 'n',
                     xlim = c(1, length(quants[, 2])),
                     ...)
      
      graphics::axis(side = 1, at = pretty_int(c(1, length(quants[, 2]))))
      
      # for (j in seq_along(x)[-1]) {
      #   graphics::lines(rowSums(pop[ , , j]),
      #                   col = 'gray')
      # }
      
      graphics::polygon(x = c(xaxs, rev(xaxs)),
                        y = c(quants[, 1], rev(quants[, 3])),
                        col = grDevices::grey(0.9),
                        border = NA)
      
      graphics::lines(quants[, 2] ~ xaxs,
                      #lwd = 2,
                      col = grDevices::grey(0.4))
      
      if (emp) {
        graphics::abline(h = round(mean(apply(pop, 3, function(x) min(rowSums(x)))), 0), lwd = 1, lty = 2)
      }
      
    }
    
    if (!is.null(stages) && stages > 0) {
      
      graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1, length(stages)))
      
      graphics::plot(pop.mn[ , stages],
                     type = 'l',
                     ylab = paste("Total Population: ", stage_names[stages]),
                     xlab = "Timesteps",
                     #lwd = 3,
                     col = graph.pal[stages],
                     xaxt = 'n',
                     xlim = c(1, length(pop.mn[ , stages])),
                     ...)
      
      graphics::axis(side = 1, at = pretty_int(c(1, length(pop.mn[ , stages]))))
      
      for (j in seq_along(x)) {
        graphics::lines(pop[ , stages, j],
                        col = 'gray',
                        lwd = 0.5)
      }
      
      graphics::lines(pop.mn[ , stages],
                      #lwd = 3,
                      col = graph.pal[stages])    
      
    }
    
  }
  
}

#' Extract objects from a 'simulation_results' object
#' 
#' The simulation results object is a list of lists containing spatial (and other) objects and
#' is organised by the following tree diagram:
#' \itemize{
#'   \item{Replicate}
#'   \itemize{
#'     \item{Timestep}
#'     \itemize{
#'       \item{Population Raster Stack}
#'       \itemize{
#'         \item{Life-Stage Raster}
#'       }
#'       \item{Habitat Suitability Raster (or Stack)}
#'       \itemize{
#'         \item{Habitat Raster (if stack is used)}
#'       }
#'       \item{Carrying Capacity Raster}
#'       \item{Other Raster Stack}
#'       \itemize{
#'         \item{Raster}
#'       }
#'       \item{...}
#'     }
#'   }
#' }
#'
#' @param x a simulation_results object 
#' @param replicate which replicate to extract from a \code{simulation_results}
#'   object
#' @param timestep timestep(s) to extract from a \code{simulation_results}
#' @param landscape_object which landscape object to extract from a
#'   \code{simulation_results} object - can be specified by name
#'   (e.g. "suitability") or index number
#' @param stage which life-stage to extract from a \code{simulation_results}
#'   object - only used for 'population' components of the landscape object
#' @param misc which misc object to extract from a \code{simulation_results}
#'   object - only used for additional spatial objects added to a landscape
#'   (e.g. disturbance layers)
#'
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)
#' 
#' pd <- population_dynamics(change = growth(egk_mat),
#'                           dispersal = kernel_dispersal(max_distance = 2000,
#'                                         dispersal_kernel = exponential_dispersal_kernel(
#'                                           distance_decay = 1000)),
#'                           density_dependence = ceiling_density())
#' 
#' sim <- simulation(landscape = ls,
#'                   population_dynamics = pd,
#'                   habitat_dynamics = NULL,
#'                   timesteps = 20)
#' 
#' # Extract the population raster for the second life-stage in the first
#' # replicate and ninth timestep
#' extract_spatial(sim, timestep = 9, stage = 2)
#' }

extract_spatial <- function (x,
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

# #' @rdname simulation
# #'
# #' @export
# #'
# #' @examples
# #'
# #' # Extract data summaries from a 'simulation_results' object
# #'
# #' extract_dispersal_info(results)
#  
# extract_dispersal_info <- function (x,
#                           replicate = 1,
#                           timestep = 1,
#                           info_object = "dispersal_failure_rate",
#                           stage = 1) {
#   om <- attr(x[[replicate]], "output_metrics")
#   
#   om[[info_object]][[timestep]][[]]
#   
# }

#' Compare minimum expected populations
#' 
#' Compare minimum expected populations from two or more 'simulation_results' objects.
#'
#' @param x a simulation_results object 
#' @param ... additional simulation results objects
#' @param interval the desired confidence interval representing the uncertainty around
#'  the expected minimum population estimates from simulation comparisons; expressed as 
#'  a whole integer between 0 and 100 (default value is 95).
#' @param all_points should the expected minimum populations from all simulation
#'  replicates be shown on the plot?
#' @param simulation_names an optional character vector of simulation names to override
#' the defaults
#'
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' ls <- landscape(population = egk_pop, suitability = egk_hab, carrying_capacity = egk_k)
#' 
#' # Create populations dynamics with and without ceiling density dependence
#' pd1 <- population_dynamics(change = growth(egk_mat),
#'                            dispersal = kernel_dispersal(max_distance = 2000,
#'                            dispersal_kernel = exponential_dispersal_kernel(distance_decay = 1000)),
#'                            density_dependence = ceiling_density())
#' pd2 <- population_dynamics(change = growth(egk_mat),
#'                            dispersal = kernel_dispersal(max_distance = 2000,
#'                            dispersal_kernel = exponential_dispersal_kernel(distance_decay = 1000)))
#' 
#' # Run first simulation with ceiling density dependence and three replicates
#' sim1 <- simulation(landscape = ls,
#'                    population_dynamics = pd1,
#'                    habitat_dynamics = NULL,
#'                    timesteps = 20,
#'                    replicates = 3)
#'                    
#' # Run second simulation without ceiling density dependence and three replicates
#' sim2 <- simulation(landscape = ls,
#'                    population_dynamics = pd2,
#'                    habitat_dynamics = NULL,
#'                    timesteps = 20,
#'                    replicates = 3)
#' 
#' compare_emp(sim1, sim2)
#' }

compare_emp <- function (x, ..., interval = 95, all_points = FALSE, simulation_names = NULL) {
  
  # read in simulation objects to compare
  sim_objects <- list(x, ...)
  n_objects <- length(sim_objects)
  
  # get names of simulations
  sim_names <- as.character(substitute(list(x, ...)))[-1L]
  
  interval_range <- c((100 - interval) / 2, 100 - (100 - interval) / 2) / 100
  
  # initiate table of values
  df <- data.frame("name" = sim_names,
                   "emp_mean" = NA,
                   "emp_lower" = NA,
                   "emp_upper" = NA)
  
  if(is.null(simulation_names)) simulation_names <- sim_names
  
  # populate table with emp mean and error values
  for (i in seq_len(n_objects)){
    pops <- get_pop_simulation(sim_objects[[i]])
    min_total_pops <- apply(pops, 3, function(x) min(rowSums(x)))
    emp_mean <- mean(min_total_pops)
    emp_lower <- stats::quantile(min_total_pops, interval_range)[1]
    emp_upper <- stats::quantile(min_total_pops, interval_range)[2]
    df[i, -1] <- c(emp_mean, emp_lower, emp_upper)
  }
  
  graphics::par(mar=c(4, 4.5, 1.5, 1.5) + 0.1)
  
  graphics::plot(NULL,
                 xlim = c(0.5, n_objects + 0.5),
                 xaxt = "n",
                 xlab = "Simulation Name",
                 ylim = range(c(df$emp_lower, df$emp_upper)),
                 yaxt = "n",
                 ylab = "",
                 main = "",
                 lwd = 0.5)
  if (all_points == TRUE) {
    for (i in seq_len(n_objects)){
      pops <- get_pop_simulation(sim_objects[[i]])
      min_total_pops <- apply(pops, 3, function(x) min(rowSums(x)))
      graphics::points(jitter(rep(i, length(min_total_pops))),
                       min_total_pops,
                       col = "lightgrey",
                       pch = 19,
                       cex = 0.8)
    }
  }
  graphics::points(seq_len(n_objects),
                   df$emp_mean,
                   pch = 19)
  graphics::arrows(seq_len(n_objects),
                   df$emp_lower,
                   seq_len(n_objects),
                   df$emp_upper,
                   length=0.05,
                   angle=90,
                   code=3)
  graphics::axis(1, at = seq_len(n_objects), labels = simulation_names)
  graphics::axis(2, at = pretty(range(c(df$emp_lower, df$emp_upper)), 5))
  graphics::mtext(paste0("Minimum Population (", interval, "% Interval)"),
                  side = 2,
                  line = 2.5)
}


##########################
### internal functions ###
##########################
as.simulation_results <- function (simulation_results) {
  as_class(simulation_results, "simulation_results", "list")
}

simulate <- function (i, landscape, population_dynamics, habitat_dynamics, timesteps, verbose, stash, is_multisession) {
  
  # if we are running in parallel, make sure the steps stash is the one from the calling session
  if (is_multisession) {
    replace_stash(stash)
  }
  
  timesteps <- seq_len(timesteps)
  if (verbose == TRUE && inherits(future::plan(), "sequential")) pb <- utils::txtProgressBar(min = 0, max = max(timesteps), style = 3)
  
  output_landscapes <- list()
  output_metrics <- list()
  output_metrics$dispersal_failure_rate <- list()

  for (timestep in timesteps) {

    for (dynamic_function in habitat_dynamics) {
      landscape <- dynamic_function(landscape, timestep)
    }
    
    # 22.01.20 - store carrying capacity function
    is_k_function <- is.function(landscape$carrying_capacity)
    
    if (is_k_function) {
      k_fun <- landscape$carrying_capacity
    }
    
    # 22.01.20 - Added to only transform carrying capacity raster once at each timestep
    landscape$carrying_capacity <- get_carrying_capacity(landscape, timestep)
    
    landscape <- population_dynamics(landscape, timestep)
    
    landscape_out <- landscape
    
    # get names of objects within the landscape object
    landscape_names <- names(landscape_out)
    
    # loop through names of objects and check if not null and greater than length one
    for (name in landscape_names) {
      
      # 22.01.20 - # if (name == "carrying_capacity" && !is.null(landscape_out[[name]]) && is.function(landscape_out[[name]])) {
        
      # 22.01.20 - #   landscape_out$carrying_capacity <- get_carrying_capacity(landscape, timestep)
        
      # 22.01.20 - # } else {
        
        if (name != "population" && !is.function(landscape_out[[name]]) &&
            !is.null(landscape_out[[name]]) && raster::nlayers(landscape_out[[name]]) > 1) {
          
          landscape_out[[name]] <- landscape_out[[name]][[timestep]]
          
        # 22.01.20 - # }
      }
    }
    
    # if (!is.null(landscape_out$suitability) && raster::nlayers(landscape_out$suitability) > 1) {
    #   landscape_out$suitability <- landscape_out$suitability[[timestep]]
    # }
    # if (!is.null(landscape_out$carrying_capacity) && is.function(landscape_out$carrying_capacity)) {
    #   landscape_out$carrying_capacity <- get_carrying_capacity(landscape, timestep)
    # }
    
    output_landscapes[[timestep]] <- landscape_out
    
    # 22.01.20 - return carrying capacity function
    if (is_k_function) {
      landscape$carrying_capacity <- k_fun
    }

    if (!is.null(steps_stash$dispersal_stats)) {
      n_stages <- length(steps_stash$dispersal_stats)
      n_stages_disperse <- which(unlist(lapply(steps_stash$dispersal_stats, function(x) !is.null(x))) == TRUE)
      dispersal_failure_rate <- replicate(n_stages, landscape$population[[1]] * 0)
      
      for (i in n_stages_disperse) {
        dispersal_failure_rate[[i]][] <- steps_stash$dispersal_stats[[i]]
      }
      
      output_metrics$dispersal_failure_rate[[timestep]] <- dispersal_failure_rate
    }
    
    if (verbose == TRUE && inherits(future::plan(), "sequential")) utils::setTxtProgressBar(pb, timestep)
    
  }
  
  if (verbose == TRUE && inherits(future::plan(), "sequential")) close(pb)
  
  as_class(output_landscapes, "replicate", "list")
  attr(output_landscapes, "output_metrics") <- output_metrics
  
  return(output_landscapes)
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
  
  pop_array <- array(dim=c(timesteps, total_stages, sims))
  
  for(i in seq_len(sims)) {
    pop_array[, , i] <- get_pop_replicate(x[[i]])
  }
  return(pop_array)
}
