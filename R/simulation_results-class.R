#' Run an simulation
#'
#' A simulation changes state objects based on selected dynamics over a
#' specified number of timesteps.
#' 
#' 
#'
#' @rdname simulation_results
#'
#' @param state a state object - static habitat, population, and demography in a timestep
#' @param dynamics a dynamics object - modules that change habitat, population, and demography during a simulation
#' @param timesteps number of timesteps used in one simulation or to display when plotting rasters
#' @param replicates number simulations to perform
#' @param parallel should parallel processors be used for simulations (default is FALSE)
#' @param x an simulation_results object
#' @param object the state object to plot - can be 'population' (default), 'habitat_suitability' or 'carrying_capacity'
#' @param type the plot type - 'graph' (default) or 'raster'
#' @param stage life-stage to plot - must be specified for 'raster' plot types; default is NULL and all life-stages will be plotted
#' @param animate if plotting type 'raster' would you like to animate the rasters as a gif?
#' @param panels the number of columns and rows to use when plotting raster timeseries - default is 3 x 3 (e.g. c(3,3) )
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
#' r <- raster(system.file("external/test.grd", package="raster"))
#'
#' mat <- matrix(c(0.000,0.000,0.302,0.302,
#'                 0.940,0.000,0.000,0.000,
#'                 0.000,0.884,0.000,0.000,
#'                 0.000,0.000,0.793,0.793),
#'               nrow = 4, ncol = 4, byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('Stage_1','Stage_2','Stage_3','Stage_4')
#'
#' pop <- stack(replicate(4, ceiling(r * 0.2)))
#'
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"),
#'                               carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(transition_matrix = mat)
#' test_population <- build_population(pop)
#'
#' test_state <- build_state(test_habitat, test_demography, test_population)
#'
#' test_dynamics <- build_dynamics(build_habitat_dynamics(),
#'                                 build_demography_dynamics(),
#'                                 build_population_dynamics())
#'
#' results <- simulation(test_state, test_dynamics, timesteps = 10, replicates = 2)

simulation <- function(state, dynamics, timesteps, replicates=1, parallel=FALSE){

  if (parallel) future::plan(multiprocess)
  simulation_results <- future.apply::future_lapply(seq_len(replicates),
                                      FUN = simulate,
                                      state = state,
                                      dynamics = dynamics,
                                      timesteps = timesteps,
                                      future.seed = FALSE)

  future::plan("default")
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

plot.simulation_results <- function (x, object = "population", type = "graph", stage = NULL, animate = FALSE, timesteps = c(1:9), panels = c(3,3), ...){
  
  stages <- raster::nlayers(x[[1]][[1]]$population$population_raster)
  stage_names <- colnames(x[[1]][[1]]$demography$global_transition_matrix)

  graph.pal <- c("#94d1c7",
                 "#cccc2b",
                 "#bebada",
                 "#fb8072",
                 "#80b1d3",
                 "#fdb462",
                 "#b3de69",
                 "#fccde5",
                 "#969696",
                 "#bc80bd"
  )
  
  if (length(x) == 1) {
   
    if (object == "population") {
      
      pop <- get_pop_replicate(x[[1]])
      
      if (type == "graph") {
        
        if (is.null(stage)) {
          
          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,stages))
          
          for (i in seq_len(stages)) {
            
            graphics::plot(pop[ , i],
                           type='l',
                           ylab=paste("Total Population: ",stage_names[i]),
                           xlab="Time (years)",
                           lwd=3,
                           col=graph.pal[i],
                           ylim=c(pretty(floor(min(pop)))[1], pretty(ceiling(max(pop)))[2]))
            
            graphics::abline(h=raster::cellStats(x[[1]][[1]]$habitat$carrying_capacity,sum)/stages,
                             lwd=1,
                             lty=2)
            
          }
          
        }
        
        if(!is.null(stage) && stage == 0) {
          
          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
          
          graphics::plot(rowSums(pop),
                         type='l',
                         ylab="Total Population (all stages)",
                         xlab="Time (years)",
                         lwd=3,
                         col="black",
                         ylim=c(pretty(floor(min(rowSums(pop))))[1], pretty(ceiling(max(rowSums(pop))))[2]))
          
          graphics::abline(h=raster::cellStats(x[[1]][[1]]$habitat$carrying_capacity,sum),
                           lwd=1,
                           lty=2)
          
        }
        
        if (!is.null(stage) && stage > 0) {
          
          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
          
          graphics::plot(pop[, stage],
                         type='l',
                         ylab=paste("Total Population: ",stage_names[stage]),
                         xlab="Time (years)",
                         lwd=3,
                         col=graph.pal[stage],
                         ylim=c(pretty(floor(min(pop[, stage])))[1], pretty(ceiling(max(pop[, stage])))[2]))
          
          graphics::abline(h=raster::cellStats(x[[1]][[1]]$habitat$carrying_capacity,sum)/stages,
                           lwd=1,
                           lty=2)
          
        }
        
      }
      
      if (type == "raster") {
        
        if (is.null(stage)) {
          stop("Please provide a life-stage when plotting \npopulation rasters or specify zero (0) for a sum of all life-stages")
        }
        
        if (length(timesteps) > 20) {
          cat ("Note, you have specified to plot rasters for more than 20 timesteps;",
               "\nthis could take a while depending on the resolution of your landscape.",
               "\nPress [enter] to continue, or [esc] to cancel and adjust settings.")
          line <- readline()
        }
        
        if(stage == 0) {
          
          rasters_sum <- raster::stack(lapply(x[[1]], function (state) sum(state$population$population_raster)))
          
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
            
            # ts <- seq_len(raster::nlayers(rasters_sum))
            # groups <- split(ts, ceiling(seq_along(ts)/9))
            # 
            # for (i in seq_along(groups)) {
            #   
            #   group <- groups[[i]]
            #   print(rasterVis::levelplot(rasters_sum[[group]],
            #                              scales = list(draw = FALSE),
            #                              margin = list(draw = FALSE),
            #                              at = breaks,
            #                              col.regions = viridisLite::viridis(length(breaks)-1),
            #                              colorkey = list(space = "bottom",
            #                                              title = "individuals",
            #                                              width = 0.6),
            #                              main="population"))
            # }
           
            print(rasterVis::levelplot(rasters_sum[[timesteps]],
                                       scales = list(draw = FALSE),
                                       margin = list(draw = FALSE),
                                       at = breaks,
                                       col.regions = viridisLite::viridis(length(breaks)-1),
                                       colorkey = list(space = "bottom",
                                                       title = "individuals",
                                                       width = 0.6),
                                       main = "population",
                                       layout = panels))
            
          }

        } else {
          
          rasters <- raster::stack(lapply(x[[1]], function (state) state$population$population_raster[[stage]]))
          
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
            
            # ts <- seq_len(raster::nlayers(rasters))
            # groups <- split(ts, ceiling(seq_along(ts)/9))
            # 
            # for (i in seq_along(groups)) {
            #   
            #   group <- groups[[i]]
            #   print(rasterVis::levelplot(rasters[[group]],
            #                              scales = list(draw = FALSE),
            #                              margin = list(draw = FALSE),
            #                              at = breaks,
            #                              col.regions = viridisLite::viridis(length(breaks)-1),
            #                              colorkey = list(space = "bottom",
            #                                              title = "individuals",
            #                                              width = 0.6),
            #                              main="population"))
            # }
            
            print(rasterVis::levelplot(rasters[[timesteps]],
                                       scales = list(draw = FALSE),
                                       margin = list(draw = FALSE),
                                       at = breaks,
                                       col.regions = viridisLite::viridis(length(breaks)-1),
                                       colorkey = list(space = "bottom",
                                                       title = "individuals",
                                                       width = 0.6),
                                       main = "population",
                                       layout = panels))
            
          }
        }
      }
    }
    
    if (object == "habitat_suitability") {
      
      rasters <- raster::stack(lapply(x[[1]], function (state) state$habitat$habitat_suitability))
      
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
        
      # ts <- seq_len(raster::nlayers(rasters))
      # groups <- split(ts, ceiling(seq_along(ts)/9))
      # 
      # for (i in seq_along(groups)) {
      #   
      #   group <- groups[[i]]
      #   print(rasterVis::levelplot(rasters[[group]],
      #                              scales = list(draw = FALSE),
      #                              margin = list(draw = FALSE),
      #                              at = breaks,
      #                              col.regions = viridisLite::viridis(length(breaks)-1),
      #                              colorkey = list(space = "bottom",
      #                                              title = "index",
      #                                              width = 0.6),
      #                              main="habitat"))
      # }
      
        print(rasterVis::levelplot(rasters[[timesteps]],
                                   scales = list(draw = FALSE),
                                   margin = list(draw = FALSE),
                                   at = breaks,
                                   col.regions = viridisLite::viridis(length(breaks)-1),
                                   colorkey = list(space = "bottom",
                                                   title = "index",
                                                   width = 0.6),
                                   main = "habitat",
                                   layout = panels))
      }
        
    }
    
    if (object == "carrying_capacity") {
      
      rasters <- raster::stack(lapply(x[[1]], function (state) state$habitat$carrying_capacity))
      
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
      
      # ts <- seq_len(raster::nlayers(rasters))
      # groups <- split(ts, ceiling(seq_along(ts)/9))
      # 
      # for (i in seq_along(groups)) {
      #   
      #   group <- groups[[i]]
      #   print(rasterVis::levelplot(rasters[[group]],
      #                              scales = list(draw = FALSE),
      #                              margin = list(draw = FALSE),
      #                              at = breaks,
      #                              col.regions = viridisLite::viridis(length(breaks)-1),
      #                              colorkey = list(space = "bottom",
      #                                              title = "individuals",
      #                                              width = 0.6),
      #                              main="k"))
      # }
      
        print(rasterVis::levelplot(rasters[[timesteps]],
                                   scales = list(draw = FALSE),
                                   margin = list(draw = FALSE),
                                   at = breaks,
                                   col.regions = viridisLite::viridis(length(breaks)-1),
                                   colorkey = list(space = "bottom",
                                                   title = "individuals",
                                                   width = 0.6),
                                   main = "k",
                                   layout = panels))
      }
        
    }
     
  }
  
  if (length(x) > 1) {
    
    pop <- get_pop_simulation(x)
    pop.mn <- round(apply(pop, c(1,2), mean), 0)
    
    if (type == "raster" | object == "habitat_suitability" | object == "carrying_capacity") {
      stop("Raster plotting is only available for single simulations")
    }
    
    if (is.null(stage)) {
      
      graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,stages))
      
      for (i in seq_len(stages)) {
        
        graphics::plot(pop.mn[, i],
                       type = 'l',
                       ylab = paste("Total Population: ",stage_names[i]),
                       xlab = "Time (years)",
                       lwd = 3,
                       col = graph.pal[i])

        for (j in seq_along(x)) {
          graphics::lines(pop[ , i, j],
                          col = 'gray')
        }
                
        #graphics::abline(h=raster::cellStats(x[[1]][[1]]$habitat$carrying_capacity,sum)/stages,
                         #lwd=1,
                         #lty=2)
        
      }
      
    }
    
    if(!is.null(stage) && stage == 0) {
      
      graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
      
      # draw the 95% CI polygon (if available) and median line
      quants <- t(apply(apply(pop, 3, rowSums),1, stats::quantile, c(0.025, 0.5, 0.975)))
      
      xaxs <- seq_len(nrow(pop[,,1]))
      
      graphics::plot(quants[,2],#rowSums(pop[ , , 1]),
                     type = 'l',
                     ylab = "Total Population (all stages)",
                     xlab = "Time (years)",
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
      
      
      graphics::abline(h=raster::cellStats(x[[1]][[1]]$habitat$carrying_capacity,sum),
                       lwd=1,
                       lty=2)
      
    }
    
    if (!is.null(stage) && stage > 0) {
      
      graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
      
      graphics::plot(pop.mn[ , stage],
                     type = 'l',
                     ylab = paste("Total Population: ",stage_names[stage]),
                     xlab = "Time (years)",
                     lwd = 3,
                     col = graph.pal[stage])
      
      for (j in seq_along(x)) {
        graphics::lines(pop[ , stage, j],
                        col = 'gray')
      }
      
      graphics::abline(h=raster::cellStats(x[[1]][[1]]$habitat$carrying_capacity,sum)/stages,
                       lwd=1,
                       lty=2)
      
    }
  
  }
  
}

##########################
### internal functions ###
##########################

as.simulation_results <- function (simulation_results) {
  as_class(simulation_results, "simulation_results", "list")
}

simulate <- function (i, state, dynamics, timesteps = 100) {
  timesteps <- seq_len(timesteps)
  output_states <- iterate_system(state, dynamics, timesteps)
  as_class(output_states, "replicate", "list")
}

iterate_system <- function (state, dynamics, timesteps) {

  output_states <- list()

  pb <- utils::txtProgressBar(min = 0, max = max(timesteps), style = 3)
  for (timestep in timesteps) {
    for (dynamic_function in dynamics) {
      state <- dynamic_function(state, timestep)
    }
    output_states[[timestep]] <- state
    utils::setTxtProgressBar(pb, timestep)
  }
  close(pb)

  output_states

}

# extract populations from a simulation
get_pop_replicate <- function(x, ...) {
  stages <- raster::nlayers(x[[1]]$population$population_raster)
  idx <- which(!is.na(raster::getValues(x[[1]]$population$population_raster[[1]])))
  pops <- lapply(x, function(x) raster::extract(x$population$population_raster, idx))
  pop_sums <- lapply(pops, function(x) colSums(x))
  pop_matrix <- matrix(unlist(pop_sums), ncol = stages, byrow = TRUE)
  return(pop_matrix)
}

get_pop_simulation <- function(x, ...) {
  stages <- raster::nlayers(x[[1]][[1]]$population$population_raster)
  timesteps <- length(x[[1]])
  sims <- length(x)
  
  pop_array <- array(dim=c(timesteps,stages,sims))
  
  for(i in seq_len(sims)) {
    pop_array[, , i] <- get_pop_replicate(x[[i]])
  }
  return(pop_array)
}
