#' Run an experiment to make spatially-explicit population projections
#'
#' @description A habitat object is used to store spatially-explicit information on habitat suitability and the carrying_capacity of a landscape.
#' It is a sub-component of a \link[dhmpr]{state} object and is modified in each timestep of an experiment.
#' 
#' @rdname experiment_results
#'
#' @param state a state object - static habitat, population, and demography in a timestep
#' @param dynamics a dynamics object - modules that change habitat, population, and demography during and experiment
#' @param timesteps number of timesteps used in the experiment
#' @param results an experiment_reults object
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
#' library(dhmpr)
#' library(raster)
#' 
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"), carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(transition_matrix = fake_transition_matrix(4), dispersal_parameters = rlnorm(1))
#' test_population <- build_population(stack(replicate(4, test_habitat$carrying_capacity * 0.2)))
#' test_state <- build_state(test_habitat, test_demography, test_population)
#' fast_approximation <- build_dynamics(no_habitat_dynamics, no_demographic_dynamics, fast_population_dynamics)
#' results <- experiment(test_state, fast_approximation, timesteps = 10)

experiment <- function (state, dynamics, timesteps = 100) {
  # check stuff
  timesteps <- seq_len(timesteps)
  output_states <- iterate_system(state, dynamics, timesteps)
  set_class(output_states, "experiment_results")
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

is.experiment_results <- function (results) {
  inherits(results, 'experiment_results')
}

#' @rdname experiment_results
#' 
#' @export
#'
#' @examples
#' 
#' print(results)

print.experiment_results <- function (results) {
  cat("This is an experiment results object, for", length(results), "timesteps")
}

#' @rdname experiment_results
#'
#' @export
#'
#' @examples
#' 
#' plot(results)

plot.experiment_results <- function (results, object = "population", type = "graph", stage = NULL, ...){

  stages <- nlayers(results[[1]]$population$population_raster)
  
  pal <- colorRampPalette(
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
  
    if (object == "population") {

      if (type == "graph") {
        
        idx <- which(!is.na(getValues(results[[1]]$population$population_raster[[1]])))
        pops <- lapply(results, function(x) extract(x$population$population_raster, idx))
        pop_sums <- lapply(pops, function(x) colSums(x))
        
        stage_names <- unlist(dimnames(results[[1]]$demography$global_transition_matrix)[1])
        
        colours <- c("#94d1c7", "#cccc2b", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#969696", "#bc80bd")

        if (is.null(stage)) {
          
          par(mfrow=c(1,stages))
          
          for (i in 1:stages) {
            
            plot(unlist(lapply(pop_sums, function(x) x[[i]])),
                 type='l',
                 ylab=paste("Total Population: ",stage_names[i]),
                 xlab="Time (years)",
                 lwd=2,
                 col=colours[i]
            )
            abline(h=cellStats(results[[1]]$habitat$carrying_capacity,sum)/stages,
                  lwd=1,
                  lty=2)
            
          }
          
<<<<<<< HEAD
        }
        
        if (!is.null(stage)) {
          if(stage==0){
=======
        } else if (stage == 0) {
>>>>>>> c54c189d5ed6fc91e1f524ee23a8a2de616cb466
          
          par(mfrow=c(1,1))

            plot(unlist(lapply(pop_sums, function(x) sum(x))),
                 type='l',
                 ylab="Total Population (all stages)",
                 xlab="Time (years)",
                 lwd=2,
                 col="black"
            )
            abline(h=cellStats(results[[1]]$habitat$carrying_capacity,sum),
                   lwd=1,
                   lty=2)
<<<<<<< HEAD
          } else{
=======

        } else {
>>>>>>> c54c189d5ed6fc91e1f524ee23a8a2de616cb466
          
          par(mfrow=c(1,1))
          
          plot(unlist(lapply(pop_sums, function(x) x[[stage]])),
               type='l',
               ylab=paste("Total Population: ",stage_names[stage]),
               xlab="Time (years)",
               lwd=2,
               col=colours[stage]
          )
          abline(h=cellStats(results[[1]]$habitat$carrying_capacity,sum)/stages,
                 lwd=1,
                 lty=2)
          
          }
        }
      }  
      
      if (type == "raster") {
        
        if (is.null(stage)) stop("Please provide a life-stage when plotting population rasters")
        
        rasters <- stack(lapply(results, function (state) state$population$population_raster[[stage]]))
        
        # Find maximum and minimum population value in raster cells for all timesteps for life-stage
        scale_max <- ceiling(max(cellStats(rasters, max)))
        scale_min <- floor(min(cellStats(rasters, min)))
        
        # Produce scale of values
        breaks <- seq(scale_min, scale_max, scale_max/10)
        
        ts <- seq_len(nlayers(rasters))
        groups <- split(ts, ceiling(seq_along(ts)/9))
        
        for (i in 1:length(groups)) {
          
          par(mar = c(0, 0, 0, 0), mfrow = c(3,3))
          group <- groups[[i]]
          plot(rasters[[group]], breaks = breaks, col = pal(length(breaks)), axes = FALSE, box = FALSE)
        
        }
      }
    
    }
  
    if (object == "habitat_suitability") {
    
      rasters <- stack(lapply(results, function (state) state$habitat$habitat_suitability))
      
      ts <- seq_len(nlayers(rasters))
      groups <- split(ts, ceiling(seq_along(ts)/9))
      
      for (i in 1:length(groups)) {
        
        par(mar = c(0, 0, 0, 0), mfrow = c(3,3))
        group <- groups[[i]]
        plot(rasters[[group]], axes = FALSE, box = FALSE)
        
      }
    
    }
  
    if (object == "carrying_capacity") {
    
      rasters <- stack(lapply(results, function (state) state$habitat$carrying_capacity))
      
      ts <- seq_len(nlayers(rasters))
      groups <- split(ts, ceiling(seq_along(ts)/9))
      
      for (i in 1:length(groups)) {
        
        par(mar = c(0, 0, 0, 0), mfrow = c(3,3))
        group <- groups[[i]]
        plot(rasters[[group]], axes = FALSE, box = FALSE)
        
      }
    
    }
  
}


##########################
### internal functions ###
##########################

iterate_system <- function (state, dynamics, timesteps) {
  
  output_states <- list()
  
  pb <- txtProgressBar(min = 0, max = max(timesteps), style = 3)
  for (timestep in timesteps) {
    for (dynamic_function in dynamics) {
      state <- dynamic_function(state, timestep)
    }
    output_states[[timestep]] <- state
    setTxtProgressBar(pb, timestep)
  }
  close(pb)
  
  output_states
  
}
