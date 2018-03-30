#' Run an experiment to make spatially-explicit population projections
#'
#' @description A habitat object is used to store spatially-explicit information on habitat suitability and the carrying_capacity of a landscape.
#' It is a sub-component of a \link[steps]{state} object and is modified in each timestep of an experiment.
#'
#' @rdname experiment_results
#'
#' @param state a state object - static habitat, population, and demography in a timestep
#' @param dynamics a dynamics object - modules that change habitat, population, and demography during and experiment
#' @param timesteps number of timesteps used in the experiment
#' @param simulations number of times to simulate the experiment
#' @param x an experiment_results or simulation object
#' @param object the state object to plot - can be 'population' (default), 'habitat_suitability' or 'carrying_capacity'
#' @param type the plot type - 'graph' (default) or 'raster'
#' @param stage life-stage to plot - must be specified for 'raster' plot types; default is NULL and all life-stages will be plotted
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{experiment_results}
#'
#' @export
#'
#' @examples
#'
#' library(steps)
#' library(raster)
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
#' test_demography <- build_demography(transition_matrix = mat,
#'                                     dispersal_parameters = rlnorm(1))
#' test_population <- build_population(pop)
#'
#' test_state <- build_state(test_habitat, test_demography, test_population)
#'
#' simple_approximation <- build_dynamics(no_habitat_dynamics,
#'                                        no_demography_dynamics,
#'                                        fast_population_dynamics)
#'
#' results <- experiment(test_state, simple_approximation, timesteps = 10)

experiment <- function (state, dynamics, timesteps = 100) {
  timesteps <- seq_len(timesteps)
  output_states <- iterate_system(state, dynamics, timesteps)
  set_class(output_states, "experiment_results")
}


#' @export
#' @noRd
`[.experiment_results` <- function(x, ..., drop=TRUE) {
  structure(NextMethod(), class=class(x))
}


#' @rdname experiment_results
#'
#' @export
#'
#' @examples
#'
#' print(results)

print.experiment_results <- function (x, ...) {
  cat("This is an experiment results object, for", length(x), "timesteps")
}


#' @rdname experiment_results
#'
#' @export
#'
#' @examples
#'
#' # Test if object is of the type 'experiment results'
#'
#' is.experiment_results(results)

is.experiment_results <- function (x) {
  inherits(x, 'experiment_results')
}

#' @rdname experiment_results
#'
#' @export
#'
#' @examples
#'
#' print(results)

print.experiment_results <- function (x, ...) {
  cat("This is an experiment results object, for", length(x), "timesteps")
}

#' @rdname experiment_results
#'
#' @export
#'
#' @examples
#'
#' plot(results)

plot.experiment_results <- function (x, object = "population", type = "graph", stage = NULL, ...){

  stages <- raster::nlayers(x[[1]]$population$population_raster)
  
  stage_names <- colnames(x[[1]]$demography$global_transition_matrix)

  ras.pal <- grDevices::colorRampPalette(
    c(
      '#440154', # dark purple
      '#472c7a', # purple
      '#3b518b', # blue
      '#2c718e', # blue
      '#21908d', # blue-green
      '#27ad81', # green
      '#5cc863', # green
      '#aadc32', # lime green
      '#fde725' # yellow
    )
  )

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

    if (object == "population") {

      pop <- get_pop_experiment(x)
      
      if (type == "graph") {

        if (is.null(stage)) {

          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,stages))
          
          for (i in seq_len(stages)) {

            graphics::plot(pop[ , i],
                 type='l',
                 ylab=paste("Total Population: ",stage_names[i]),
                 xlab="Time (years)",
                 lwd=2,
                 col=graph.pal[i]
            )
            graphics::abline(h=raster::cellStats(x[[1]]$habitat$carrying_capacity,sum)/stages,
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
               lwd=2,
               col="black"
          )
          graphics::abline(h=raster::cellStats(x[[1]]$habitat$carrying_capacity,sum),
                 lwd=1,
                 lty=2)

        }

        if (!is.null(stage) && stage > 0) {

          graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))

          graphics::plot(pop[, stage],
               type='l',
               ylab=paste("Total Population: ",stage_names[stage]),
               xlab="Time (years)",
               lwd=2,
               col=graph.pal[stage]
          )
          graphics::abline(h=raster::cellStats(x[[1]]$habitat$carrying_capacity,sum)/stages,
                 lwd=1,
                 lty=2)

        }

      }

      if (type == "raster") {

        if (is.null(stage)) stop("Please provide a life-stage when plotting population rasters")

        rasters <- raster::stack(lapply(x, function (state) state$population$population_raster[[stage]]))

        # Find maximum and minimum population value in raster cells for all timesteps for life-stage
        scale_max <- ceiling(max(raster::cellStats(rasters, max)))
        scale_min <- floor(min(raster::cellStats(rasters, min)))

        # Produce scale of values
        breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)

        ts <- seq_len(raster::nlayers(rasters))
        groups <- split(ts, ceiling(seq_along(ts)/9))

        for (i in seq_along(groups)) {

          group <- groups[[i]]
          print(rasterVis::levelplot(rasters[[group]],
                          scales = list(draw = FALSE),
                          at = breaks,
                          col.regions = ras.pal(length(breaks)-1),
                          colorkey = list(space = "bottom",
                                          title = "individuals",
                                          width = 0.6
                                          ),
                          main="population"
                          )
                )
        }

      }

    }

    if (object == "habitat_suitability") {

      rasters <- raster::stack(lapply(x, function (state) state$habitat$habitat_suitability))

      # Find maximum and minimum population value in raster cells for all timesteps for life-stage
      scale_max <- ceiling(max(raster::cellStats(rasters, max)))
      scale_min <- floor(min(raster::cellStats(rasters, min)))

      # Produce scale of values
      breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)

      ts <- seq_len(raster::nlayers(rasters))
      groups <- split(ts, ceiling(seq_along(ts)/9))

      for (i in seq_along(groups)) {

        group <- groups[[i]]
        print(rasterVis::levelplot(rasters[[group]],
                        scales = list(draw = FALSE),
                        at = breaks,
                        col.regions = ras.pal(length(breaks)-1),
                        colorkey = list(space = "bottom",
                                        title = "index",
                                        width = 0.6
                                        ),
                        main="habitat"
                        )
              )
      }

    }

    if (object == "carrying_capacity") {

      rasters <- raster::stack(lapply(x, function (state) state$habitat$carrying_capacity))

      # Find maximum and minimum population value in raster cells for all timesteps for life-stage
      scale_max <- ceiling(max(raster::cellStats(rasters, max)))
      scale_min <- floor(min(raster::cellStats(rasters, min)))

      # Produce scale of values
      breaks <- seq(scale_min, scale_max, (scale_max-scale_min)/100)

      ts <- seq_len(raster::nlayers(rasters))
      groups <- split(ts, ceiling(seq_along(ts)/9))

      for (i in seq_along(groups)) {

        group <- groups[[i]]
        print(rasterVis::levelplot(rasters[[group]],
                        scales = list(draw = FALSE),
                        at = breaks,
                        col.regions = ras.pal(length(breaks)-1),
                        colorkey = list(space = "bottom",
                                        title = "individuals",
                                        width = 0.6
                                        ),
                        main="k"
                        )
              )
      }

    }

}

#' @rdname experiment_results
#'
#' @name simulation
#'
#' @export
#'
#' @importFrom future plan multiprocess future values
#'
#' @return An object of class \code{simulation_results} which
#' contains n \code{experiment_results}
#'
#' @examples
#'
#' library(future)
#' plan(multiprocess)
#' sim_results <- simulation(test_state, simple_approximation,
#'                           timesteps = 10, simulations=10)

simulation <- function(state, dynamics, timesteps, simulations){

  simulation_results <- list()
  for(ii in seq_len(simulations)){
    simulation_results[[ii]] <- future({
      experiment(state,dynamics,timesteps)
    },
    globals = list(state = state,
                   dynamics = dynamics,
                   timesteps = timesteps,
                   experiment = steps::experiment)
    )
  }
  set_class(future::values(simulation_results), "simulation_results")
}

#' @rdname experiment_results
#'
#' @export
#'
#' @examples
#'
#' plot(sim_results)

plot.simulation_results <- function (x, stage = NULL, ...){
  
  stages <- raster::nlayers(x[[1]][[1]]$population$population_raster)
  
  stage_names <- colnames(x[[1]]$demography$global_transition_matrix)
  
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
  
  pop <- get_pop_simulation(x)
  
  if (is.null(stage)) {
    
    graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,stages))
    
    for (i in seq_len(stages)) {
      graphics::plot(pop[ , i, 1],
                     type = 'l',
                     ylab = paste("Total Population: ",stage_names[i]),
                     xlab = "Time (years)",
                     lwd = 2,
                     col = graph.pal[i]
      )
      
      for (j in seq_along(x)[-1]) {
        graphics::lines(pop[ , i, j],
                        col = 'gray'
        )
      }
      
      graphics::abline(h=raster::cellStats(x[[1]][[1]]$habitat$carrying_capacity,sum)/stages,
                       lwd=1,
                       lty=2)
      
    }
    
  }
  
  if(!is.null(stage) && stage == 0) {
    
    graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
    
    graphics::plot(rowSums(pop[ , , 1]),
                   type = 'l',
                   ylab = paste("Total Population: ",stage_names[i]),
                   xlab = "Time (years)",
                   lwd = 2,
                   col = 'black'
    )
    
    for (j in seq_along(x)[-1]) {
      graphics::lines(rowSums(pop[ , , j]),
                      col = 'black'
      )
    }
    
    graphics::abline(h=raster::cellStats(x[[1]][[1]]$habitat$carrying_capacity,sum),
                     lwd=1,
                     lty=2)
    
  }
  
  if (!is.null(stage) && stage > 0) {
    
    graphics::par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
    
    graphics::plot(pop[ , stage, 1],
                   type = 'l',
                   ylab = paste("Total Population: ",stage_names[i]),
                   xlab = "Time (years)",
                   lwd = 2,
                   col = graph.pal[stage]
    )
    
    for (j in seq_along(x)[-1]) {
      graphics::lines(pop[ , stage, j],
                      col = 'gray'
      )
    }
    
    graphics::abline(h=raster::cellStats(x[[1]][[1]]$habitat$carrying_capacity,sum)/stages,
                     lwd=1,
                     lty=2)
    
  }
  
}

##########################
### internal functions ###
##########################

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

# extract populations from an experiment
get_pop_experiment <- function(x, ...) {
  idx <- which(!is.na(raster::getValues(x[[1]]$population$population_raster[[1]])))
  pops <- lapply(x, function(x) raster::extract(x$population$population_raster, idx))
  pop_sums <- lapply(pops, function(x) colSums(x))
  pop_matrix <- matrix(unlist(pop_sums), ncol = 4, byrow = TRUE)
  return(pop_matrix)
}

get_pop_simulation <- function(x, ...) {
  stages <- raster::nlayers(x[[1]][[1]]$population$population_raster)
  timesteps <- length(x[[1]])
  sims <- length(x)
  
  pop_array <- array(dim=c(timesteps,stages,sims))
  
  for(i in seq_len(sims)) {
    pop_array[, , i] <- get_pop_experiment(x[[i]])
  }
  
  return(pop_array)
}
