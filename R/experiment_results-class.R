#' Run an experiment to make spatially-explicit population projections
#'
#' @param state A state object - static habitat, population, and demography in a timestep
#' @param dynamics A dynamics object - modules that change habitat, population, and demography during and experiment
#' @param timesteps Number of timesteps used in the experiment
#' 
#' @return An object of class \code{experiment_results}
#' @export
#'
#' @examples
#' 
#' library(raster)
#' library(dhmpr)
#' 
#' r <- raster(system.file("external/test.grd", package="raster"))
#' 
#' test_habitat <- build_habitat(habitat_suitability = r / cellStats(r, "max"), carrying_capacity = ceiling(r * 0.1))
#' test_demography <- build_demography(fake_transition_matrix(4), rlnorm(1))
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


#' Print details of a experiment object
#'
#' @param results an object to print or test as an experiment_results object
#' @param ... further arguments passed to or from other methods
#'
#' @export
#'
# @examples
# results <- experiment(test_state, fast_approximation, timesteps = 10)
# print(results)

print.experiment_results <- function (results) {
  #cat("This is an experiment results object, for", length(x), "timesteps")
}

#' Plot an experiment object
#'
#' @param results 
#' @param object 
#' @param type 
#' @param ... 
#'
#' @export
#'
# @examples
# results <- experiment(test_state, fast_approximation, timesteps = 10)
# plot(results)

plot.experiment_results <- function (results, object = "population", type = "raster", stage = NULL, ...) {

  stages <- nlayers(results[[1]]$population$population_raster)
  
    if (object == "population") {

      if (type == "graph") {
        
        idx <- which(!is.na(getValues(results[[1]]$population$population_raster[[1]])))
        pops <- lapply(results, function(x) extract(x$population$population_raster, idx))
        pop_sums <- lapply(pops, function(x) colSums(x))
        
        stage_names <- unlist(dimnames(results[[1]]$demography$transition_matrix)[1])
        
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
          
        }else{
          
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
      
      if (type == "raster") {
        
        if (is.null(stage)) stop("Please provide a life-stage when plotting population rasters")
        
        rasters <- stack(lapply(results, function (state) state$population$population_raster[[stage]]))
        
        ts <- seq_len(nlayers(rasters))
        groups <- split(ts, ceiling(seq_along(ts)/9))
        
        for (i in 1:length(groups)) {
          
          par(mfrow=c(3,3))
          group <- groups[[i]]
          plot(rasters[[group]])
        
        }
      }
    
    }
  
    if (object == "habitat_suitability") {
    
      rasters <- stack(lapply(results, function (state) state$habitat$habitat_suitability))
      
      ts <- seq_len(nlayers(rasters))
      groups <- split(ts, ceiling(seq_along(ts)/9))
      
      for (i in 1:length(groups)) {
        
        par(mfrow=c(3,3))
        group <- groups[[i]]
        plot(rasters[[group]])
        
      }
    
    }
  
    if (object == "carrying_capacity") {
    
      rasters <- stack(lapply(results, function (state) state$habitat$carrying_capacity))
      
      ts <- seq_len(nlayers(rasters))
      groups <- split(ts, ceiling(seq_along(ts)/9))
      
      for (i in 1:length(groups)) {
        
        par(mfrow=c(3,3))
        group <- groups[[i]]
        plot(rasters[[group]])
        
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
