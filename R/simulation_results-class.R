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
#' @importFrom future plan multisession future value
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
    future.apply::future_lapply(..., future.globals = future.globals, future.seed = TRUE)
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
  
  # store intitial population, habitat, and carrying capacity objects
  initial_population <- landscape$population
  initial_suitability <- landscape$suitability
  initial_carrying_capacity <- get_carrying_capacity(landscape, 1)
  
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
  
  # add the initial population, habitat, and carrying capacity objects as attributes
  attr(simulation_results, "initial_population") <- initial_population
  attr(simulation_results, "initial_suitability") <- initial_suitability
  attr(simulation_results, "initial_carrying_capacity") <- initial_carrying_capacity
  
  as.simulation_results(simulation_results)
}

#' @export
#' @noRd
`[.simulation_results` <- function(x, ..., drop = TRUE) {
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
#' and spatial-explicit grids are generated for all timesteps to illustrate
#' population changes through time and space. Note, this function can be wrapped
#' in a *png()* call to write several images to disk for creating animations.
#' 
#' @param x a simulation_results object
#' @param replicates which replicates to plot (default is one, or the first)
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
#' # Plot the spatial distributions of total cell populations
#' plot(sim) 
#' }

plot.simulation_results <- function (x,
                                     replicates = 1,
                                     ...){
  
  # avoid a persistent effect on the graphics device
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  
  pop_data_stages <- get_pop_simulation(x)
  pop_data_totals <- apply(pop_data_stages, 3, function(x) rowSums(x))
  
  timesteps <- seq_len(length(x[[1]]))
  
  pop_spatial <- lapply(replicates,
                        function(r) lapply(timesteps,
                                           function(t) sum(x[[r]][[t]][[1]])))
  max_ind <- max(sapply(replicates,
                        function(r) sapply(timesteps,
                                           function(t) max(pop_spatial[[r]][[t]][],
                                                           na.rm = TRUE))))
  
  # rounded_max_ind <- 50 * ceiling(max_ind / 50)
  # cuts <- seq(0, rounded_max_ind, 50)
  
  if (max_ind < 5) {
    cuts <- pretty(seq_len(max_ind), 2)
  }else{
    cuts <- pretty(seq_len(max_ind), 5)
  }
  # cols <- c("lightgray", viridis(length(cuts) - 1))

  cols <- c("lightgray", viridis(max_ind - 1))
  layer_names <- paste0("Population in year ", timesteps)
  
  for (i in replicates){
    for (j in timesteps) {

      graphics::layout(matrix(c(1, 1, 1, 2, 2, 3), ncol = 2), width = c(1, 1, 1), height = c(3, 2, 1))
      graphics::par(mar = c(4, 3, 3, 0), mgp = c(2, 0.7, 0))
      
      graphics::plot(x = 0:j, y = c(pop_data_totals[1, i], pop_data_totals[-1, i][1:j]),
                     type = 'l',
                     ylab  = 'Population size',
                     ylim = c(0, max(pop_data_totals) * 1.1),
                     xlab = "Year",
                     xlim = c(0, max(timesteps)),
                     main = paste0("Simulation replicate ", i),
                     font.main = 1,
                     cex.main = 0.9)
      graphics::grid()
      
      mat <- t(apply(raster::as.matrix(pop_spatial[[i]][[j]]), 2, rev))
      
      graphics::par(mar = c(2, 2, 1.2, 1.2))
      graphics::image(mat,
                      col = cols,
                      axes = FALSE,
                      asp = 1,
                      zlim = c(0, max_ind))
      graphics::title(paste0(layer_names[j]),
                      line = -0.7,
                      font.main = 1,
                      cex.main = 0.9)

      legend_image <- raster::as.raster(matrix(cols, nrow = 1))
      
      graphics::plot.new()
      graphics::par(mar = c(2, 2, 1.2, 0.8))
      graphics::rasterImage(legend_image, 0, 0.8, 0.9, 1)
      graphics::text(x = (cuts / max_ind) * 0.9, y = 0.4, labels = cuts)
    }
  }
  
}


#' Extract spatial object from a 'simulation_results' object
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
#' @param timestep which timestep to extract from a \code{simulation_results}
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
#' extract_spatial(sim, replicate = 1, timestep = 9, stage = 2)
#' }

extract_spatial <- function (x,
                             replicate = 1,
                             timestep = 1,
                             landscape_object = "population",
                             stage = 1,
                             misc = 1) {
  if (landscape_object == 1 | landscape_object == "population") {
    return(x[[replicate]][[timestep]][[landscape_object]][[stage]])
  }
  if (landscape_object > 3 |
      !landscape_object == "population" |
      !landscape_object == "suitability" |
      !landscape_object == "carrying_capacity") {
    return(x[[replicate]][[timestep]][[landscape_object]][[misc]])
  }
  return(x[[replicate]][[timestep]][[landscape_object]])
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
