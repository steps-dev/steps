#' Visualise the results of a *steps* simulation
#'
#' Visualising the results of a simulation is important to verify parameter assumptions
#' and quantitative model behaviour. Both linear graphs indicating trends and
#' spatial-explicit grids containing spatial arrangement of information can be generated
#' to illustrate changes through time and space for populations, carrying capacity,
#' and habitat suitability. The expected minimum populations (EMP) can also be compared
#' for several different simulations.
#'
#' For plotting trends, see:
#' \itemize{
#' \item{\code{\link[steps]{plot_pop_trend}} to examine population changes}
#' \item{\code{\link[steps]{plot_k_trend}} to examine carrying capacity changes}
#' }
#' 
#' For plotting spatial information, see:
#' \itemize{
#' \item{\code{\link[steps]{plot_pop_spatial}} to examine population changes}
#' \item{\code{\link[steps]{plot_k_spatial}} to examine carrying capacity changes}
#' \item{\code{\link[steps]{plot_hab_spatial}} to examine habitat suitability changes}
#' }
#' 
#' For plotting and comparing expected minimum populations, see:
#' \itemize{
#' \item{\code{\link[steps]{compare_emp}} to examine how different simulations compare}
#' }
#' 
#' @name visualisation
NULL


#' Plot population trend
#' 
#' Plot linear graphs to illustrate population changes through time.
#' 
#' @param x a simulation_results object
#' @param stages life-stages to plot - by default all life-stages will be considered. 
#' Set to zero for totals (i.e. sum of all life-stages).
#' @param emp (TRUE/FALSE) add a dashed line indicating the expected minimum
#' population of the simulation (for multiple replicates only)
#' @param return_data (TRUE/FALSE) should the data used to create the plots be returned?
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
#' plot_pop_trend(sim)
#' 
#' # Plot the total population trajectory
#' plot_pop_trend(sim, stages = 0) 
#' }

plot_pop_trend <- function (x,
                            stages = NULL,
                            emp = FALSE,
                            return_data = FALSE,
                            ...){
  
  # avoid a persistent effect on the graphics device
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  
  graph_pal <- c(
    "#6da36b",
    "#eb7d75",
    "#80b1d3",
    "#bebada",
    "#f0ab7e",
    "#969696",
    "#d1a954",
    "#81b84d",
    "#ba5059",
    "#3c3c7d"
  )
  
  total_stages <- raster::nlayers(x[[1]][[1]]$population)
  stage_names <- names(x[[1]][[1]]$population)
  replicates <- length(x)
  timesteps <- length(x[[1]])
  
  pop <- get_pop_simulation(x)
  pop_mn <- round(apply(pop, c(1,2), mean), 0)
  
  if (is.null(stages)) {
    stages <- seq_len(total_stages)
  }
  
  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1, length(stages)))
  
  y_label <- paste0("Total Population: ", stage_names[stages])
  y_range <- range(pretty(pop))
  
  if (all(stages == 0)) {
    stages <- 1
    emp_0 <- TRUE
    pop <- structure(sapply(seq_len(replicates), FUN = function(x) rowSums(pop[ , , x])),
                     dim = c(timesteps + 1, 1, replicates))
    pop_mn <- data.frame(round(apply(pop, 1, mean), 0))
    emp_value <- round(mean(apply(pop, 3, function(x) min(rowSums(x)))), 0)
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1, 1))
    y_label <- "Total Population (all stages)"
    y_range <- range(pretty(pop_mn[ , 1]))
    graph_pal <- "black"
  }
  
  for (i in stages) {
    
    graphics::plot(c(0, seq_len(timesteps)),
                   pop_mn[ , i],
                   type = 'l',
                   ylab = y_label[i],
                   xlab = "Timesteps",
                   #lwd = 3,
                   col = graph_pal[i],
                   ylim = y_range,
                   xaxt = 'n',
                   xlim = c(0, timesteps),
                   ...)
    graphics::axis(side = 1, at = pretty(c(0, timesteps)))
    graphics::grid()
    
    if(replicates > 1){
      for (j in seq_len(replicates)) {
        graphics::lines(c(0, seq_len(timesteps)),
                        pop[ , i, j],
                        col = 'gray',
                        lwd = 0.5)
      }
      graphics::lines(c(0, seq_len(timesteps)),
                      pop_mn[ , i],
                      #lwd = 3,
                      col = graph_pal[i])
      if (emp & exists("emp_0")) {
        graphics::abline(h = emp_value, lwd = 1, lty = 2)
      }
    }
    
  }
  
  if (return_data == TRUE) return(pop)
}

# test <- plot_pop_trend(egk_results, emp = TRUE, return_data = TRUE)
# plot_pop_trend(egk_results, stages = c(1, 2), emp = TRUE)
# plot_pop_trend(egk_results, stages = 0, emp = TRUE)
# 
# egk_results2 <- egk_results[1]
# attr(egk_results2, "initial_population") <- attr(egk_results, "initial_population")
# 
# plot_pop_trend(egk_results2, emp = TRUE)
# plot_pop_trend(egk_results2, stages = c(1, 2), emp = TRUE)
# plot_pop_trend(egk_results2, stages = 0, emp = TRUE)


#' Plot carrying capacity (k) trend
#' 
#' Plot linear graphs to illustrate carrying capacity changes through time.
#' 
#' @param x a simulation_results object
#' @param summary_stat how to summarize the values across the landscape - "mean" (default)
#' or "sum"
#' #' @param return_data (TRUE/FALSE) should the data used to create the plots be returned?
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
#' # Plot the carrying capacity trajectories
#' plot_k_trend(sim)
#' }

plot_k_trend <- function (x,
                          summary_stat = "mean",
                          return_data = FALSE,
                          ...){
  
  # avoid a persistent effect on the graphics device
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  
  replicates <- length(x)
  timesteps <- length(x[[1]])
  
  init_k <- rep(round(mean(attr(x, "initial_carrying_capacity")[], na.rm = TRUE)), replicates)
  k <- sapply(seq_len(replicates),
              FUN = function(r) sapply(seq_len(timesteps),
                                       FUN = function(t) round(mean(x[[r]][[t]][[3]][], na.rm = TRUE))))
  
  if (summary_stat == "sum"){
    init_k <- rep(round(sum(attr(x, "initial_carrying_capacity")[], na.rm = TRUE)), replicates)
    k <- sapply(seq_len(replicates),
                FUN = function(r) sapply(seq_len(timesteps),
                                         FUN = function(t) round(sum(x[[r]][[t]][[3]][], na.rm = TRUE))))
  }
  
  all_k <- rbind(init_k, k)
  all_k_mn <- round(apply(all_k, 1, mean), 0)
  
  y_range <- range(pretty(all_k))
  
  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1, 1))
  
  graphics::plot(c(0, seq_len(timesteps)),
                 all_k_mn,
                 type = 'l',
                 ylab = paste0("Landscape Carrying Capacity (", summary_stat, ")"),
                 xlab = "Timesteps",
                 #lwd = 3,
                 ylim = y_range,
                 xaxt = 'n',
                 xlim = c(0, timesteps),
                 ...)
  graphics::axis(side = 1, at = pretty(c(0, timesteps)))
  
  if(replicates > 1){
    for (j in seq_len(replicates)) {
      graphics::lines(c(0, seq_len(timesteps)),
                      all_k[ , j],
                      col = 'gray',
                      lwd = 0.5)
    }
    graphics::lines(c(0, seq_len(timesteps)),
                    all_k_mn)
  }
  
  if (return_data == TRUE) return(all_k)
}

# test <- plot_k_trend(egk_results, return_data = TRUE)
# plot_k_trend(egk_results, summary_stat = "sum")
# 
# attr(egk_results2, "initial_carrying_capacity") <- attr(egk_results, "initial_carrying_capacity")
# 
# plot_k_trend(egk_results2)
# plot_k_trend(egk_results2, summary_stat = "sum")


#' Plot population spatial information
#' 
#' Plot spatial grids to illustrate population changes through time.
#' 
#' @param x a simulation_results object
#' @param stage life-stage to plot - defaults to totals of all life stages. 
#' Set to zero for totals (i.e. sum of all life-stages).
#' @param replicate replicate to plot - note, only one replicate can be plotted
#' at a time. The default is to plot the first replicate
#' @param timesteps timesteps to plot
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
#' plot_pop_spatial(sim)
#' }

plot_pop_spatial <- function (x,
                              stage = 0,
                              replicate = 1,
                              timesteps = NULL,
                              ...){
  
  total_stages <- raster::nlayers(x[[1]][[1]]$population)
  stage_names <- names(x[[1]][[1]]$population)
  total_timesteps <- length(x[[1]])
  
  if (is.null(timesteps)){
    timesteps <- seq_len(total_timesteps)
  }
  
  r <- lapply(seq_len(total_timesteps), FUN = function(t) raster::stack(lapply(seq_len(total_stages), FUN = function(s) extract_spatial(x,
                                                                                                                                        replicate = replicate,
                                                                                                                                        timestep = t,
                                                                                                                                        stage = s))))
  
  if(stage == 0){
    r <- raster::stack(lapply(seq_len(total_timesteps), FUN = function(t) sum(r[[t]])))
    names(r) <- paste0("total_", seq_len(total_timesteps))
  }else{
    r <- raster::stack(lapply(seq_len(total_timesteps), FUN = function(t) r[[t]][[stage]]))
    names(r) <- paste0(stage_names[stage], "_", seq_len(total_timesteps))
  }
  
  r <- r[[timesteps]]
  
  plot_spatial(raster_stack = r,
               label = "Individuals",
               ...)
}


#' Plot carrying capacity spatial information
#' 
#' Plot spatial grids to illustrate carrying capacity changes through time.
#' 
#' @param x a simulation_results object.
#' @param replicate replicate to plot - note, only one replicate can be plotted
#' at a time. The default is to plot the first replicate
#' @param timesteps timesteps to plot
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
#' plot_k_spatial(sim)
#' }

plot_k_spatial <- function (x,
                            replicate = 1,
                            timesteps = NULL,
                            ...){
  
  total_timesteps <- length(x[[1]])
  
  if (is.null(timesteps)){
    timesteps <- seq_len(total_timesteps)
  }
  
  r <- raster::stack(lapply(seq_len(total_timesteps), FUN = function(t) extract_spatial(x,
                                                                                        replicate = replicate,
                                                                                        timestep = t,
                                                                                        landscape_object = "carrying_capacity")))
  names(r) <- paste0("Timestep_", seq_len(total_timesteps))
  
  r <- r[[timesteps]]
  
  plot_spatial(raster_stack = r,
               label = "Carrying Capacity",
               ...)
}


#' Plot habitat suitability spatial information
#' 
#' Plot spatial grids to illustrate habitat suitability changes through time.
#' 
#' @param x a simulation_results object.
#' @param replicate replicate to plot - note, only one replicate can be plotted
#' at a time. The default is to plot the first replicate
#' @param timesteps timesteps to plot
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
#' plot_hab_spatial(sim)
#' }

plot_hab_spatial <- function (x,
                              replicate = 1,
                              timesteps = NULL,
                              ...){
  
  total_timesteps <- length(x[[1]])
  
  if (is.null(timesteps)){
    timesteps <- seq_len(total_timesteps)
  }
  
  r <- raster::stack(lapply(seq_len(total_timesteps), FUN = function(t) extract_spatial(x,
                                                                                        replicate = replicate,
                                                                                        timestep = t,
                                                                                        landscape_object = "suitability")))
  names(r) <- paste0("Timestep_", seq_len(total_timesteps))
  
  r <- r[[timesteps]]
  
  plot_spatial(raster_stack = r,
               label = "Habitat Suitability",
               ...)
}


#' Compare minimum expected populations
#' 
#' Compare minimum expected populations from two or more 'simulation_results' objects.
#'
#' @param x a simulation_results object 
#' @param ... additional simulation results objects
#' @param show_interval should the interval bars be shown on the plot? Default is TRUE.
#' @param interval the desired confidence interval representing the uncertainty around
#'  the expected minimum population estimates from simulation comparisons; expressed as 
#'  a whole integer between 0 and 100 (default value is 95).
#' @param all_points should the expected minimum populations from all simulation
#'  replicates be shown on the plot? Default is FALSE.
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
#'                            dispersal = kernel_dispersal(max_distance = 1000,
#'                            dispersal_kernel = exponential_dispersal_kernel(distance_decay = 500)),
#'                            density_dependence = ceiling_density())
#' pd2 <- population_dynamics(change = growth(egk_mat),
#'                            dispersal = kernel_dispersal(max_distance = 3000,
#'                            dispersal_kernel = exponential_dispersal_kernel(distance_decay = 1500)))
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

compare_emp <- function (x, ..., show_interval = TRUE, interval = 95, all_points = FALSE, simulation_names = NULL) {
  
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
  if (show_interval == TRUE) {
    graphics::arrows(seq_len(n_objects),
                     df$emp_lower,
                     seq_len(n_objects),
                     df$emp_upper,
                     length=0.05,
                     angle=90,
                     code=3)
  }
  graphics::axis(1, at = seq_len(n_objects), labels = simulation_names)
  graphics::axis(2, at = pretty(range(c(df$emp_lower, df$emp_upper)), 5))
  graphics::mtext(paste0("Minimum Population (", interval, "% Interval)"),
                  side = 2,
                  line = 2.5)
}


##########################
### internal functions ###
##########################

plot_spatial <- function(raster_stack, label, ...){
  
  # avoid a persistent effect on the graphics device
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  
  scale_max <- ceiling(max(raster::cellStats(raster_stack, max)))
  scale_min <- floor(min(raster::cellStats(raster_stack, min)))
  
  breaks <- seq(scale_min, scale_max, (scale_max - scale_min) / 100)
  
  ifelse(any(raster_stack[] == 0, na.rm = TRUE),
         colour_range <- c("#bfbfbfff", viridisLite::viridis(length(breaks) - 1)),
         colour_range <- viridisLite::viridis(length(breaks) - 1))
  
  graphics::par(mar = c(2, 0, 0, 0), mfrow = c(1,1))
  print(rasterVis::levelplot(raster_stack,
                             scales = list(draw = FALSE),
                             margin = list(draw = FALSE),
                             at = breaks,
                             col.regions = colour_range,
                             colorkey = list(space = "bottom",
                                             width = 0.4),
                             par.settings = list(layout.heights = list(xlab.key.padding = 1),
                                                 strip.background = list(col = "#ffffff")),
                             xlab = label,
                             layout = c(3, 3),
                             ...))
  
}
