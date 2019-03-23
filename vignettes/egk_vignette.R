## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(dpi=300,fig.width=7)

## ---- message = FALSE----------------------------------------------------
library(steps)
library(raster)
library(viridisLite)
library(future)

## ---- message = FALSE----------------------------------------------------

egk_mat <- matrix(c(0.00,0.00,1.00,
                    0.50,0.00,0.00,
                    0.00,0.85,0.85),
                  nrow = 3,
                  ncol = 3,
                  byrow = TRUE)
colnames(egk_mat) <- rownames(egk_mat) <- c('juvenile','subadult','adult')

egk_mat_stoch <- matrix(c(0.00,0.00,0.20,
                          0.05,0.00,0.00,
                          0.00,0.10,0.05),
                  nrow = 3,
                  ncol = 3,
                  byrow = TRUE)
colnames(egk_mat_stoch) <- rownames(egk_mat_stoch) <- c('juvenile','subadult','adult')

## ---- message = FALSE, fig.align = "center", out.width = "100%"----------
egk_hab

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(egk_hab, box = FALSE, axes = FALSE, col = viridis(100), main = "Habitat Suitability")

## ---- message = FALSE, fig.align = "center", out.width = "100%"----------
egk_k

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(egk_k, box = FALSE, axes = FALSE, col = viridis(100))

## ---- message = FALSE, fig.align = "center", out.width = "100%"----------
egk_pop

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spplot(egk_pop, col.regions = viridis(100))

## ---- message = FALSE----------------------------------------------------
library(foreach)
stable_states <- abs( eigen(egk_mat)$vectors[,1] / sum(eigen(egk_mat)$vectors[,1]))
popN <- stack(replicate(ncol(egk_mat), egk_k)) * stable_states
idx <- which(!is.na(getValues(popN[[1]])))
pop <- stack(
  foreach(i = 1:nlayers(popN)) %do% {
    max_pop <- ceiling(cellStats(popN[[i]], max, na.rm = T))
    pop_values <- popN[[i]][idx]
    popN[[i]][idx] <- rbinom(prob = (pop_values/max_pop), size = max_pop, n = length(pop_values))
    popN[[i]]
  })
names(pop) <- colnames(egk_mat)
pop

## ---- message = FALSE----------------------------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = NULL,
                           carrying_capacity = NULL)

## ---- message = FALSE----------------------------------------------------
egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                        dispersal = NULL,
                                        modification = NULL,
                                        density_dependence = NULL)

## ---- message = FALSE,  results = 'hide', progress = FALSE---------------
egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 1,
                          verbose = FALSE)

## ---- message = FALSE, fig.align = "center", out.width = "100%"----------
plot(egk_results)

## ---- message = FALSE, fig.width = 4, fig.align = "center"---------------
plot(egk_results, stage = 0)

## ---- message = FALSE, fig.width = 4, fig.align = "center"---------------
plot(egk_results, stage = 2, newplot = TRUE)

## ---- message = FALSE, fig.align = "center", out.width = "100%", fig.show='hold'----
plot(egk_results, type = "raster", stage = 2, timesteps = c(1, 10, 20), panels = c(3, 1))

## ---- message = FALSE, eval = FALSE--------------------------------------
#  plot(egk_results, type = "raster", stage = 2, timesteps = c(1, 10, 20), animate = TRUE)

## ---- message = FALSE, progress = FALSE----------------------------------
plan(multiprocess, workers = 3) # This is how we specify to simulate replicates on separate processors in parallel

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, progress = FALSE, fig.align = "center", out.width = "100%"----
egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        global_stochasticity = egk_mat_stoch),
                                        dispersal = NULL,
                                        modification = NULL,
                                        density_dependence = NULL)

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          demo_stochasticity = "none",
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, progress = FALSE, fig.align = "center", out.width = "100%"----
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = NULL)

egk_pop_dynamics <- population_dynamics(
  change = growth(transition_matrix = egk_mat,
                  transition_function = modified_transition(egk_mat,
                                                            survival_layer = "suitability",
                                                            fecundity_layer = "suitability")),
  dispersal = NULL,
  modification = NULL,
  density_dependence = NULL)

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, eval = FALSE--------------------------------------
#  deterministic_transitions <- function(transition_matrix) {
#  
#    dim <- nrow(transition_matrix)
#  
#    function (landscape, timestep) {
#  
#      #### This assumes that the files are named with a particular convention
#      #### and will not work for all cases. User must change code below accordingly.
#  
#        # get metrics and constructor info
#        cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
#        current_timestep <- sprintf("%02i", timestep)
#        n_cells <- length(which(!is.na(raster::getValues(landscape$population[[1]]))))
#  
#        # get relevant rasters - note, your working directory will be different
#        files <- list.files("../working/rasters", pattern = paste0("_", current_timestep, "_"))
#  
#        #initialise array
#        transition_array <- array(0, dim = c(dim, dim, n_cells))
#  
#        # populate array:
#        for (file in files) {
#          r <- as.integer(substr(substr(file, nchar(file) - (8-1), nchar(file)), 1, 2))
#          c <- as.integer(substr(substr(file, nchar(file) - (6-1), nchar(file)), 1, 2))
#          #note, your working directory will be different
#          transition_array[r, c, ] <- raster::raster(paste0("../working/rasters/",
#                                                            file))[cell_idx]
#      }
#  
#      #### Return array with required dimensions
#      transition_array
#    }
#  }
#  
#  egk_pop_dynamics <- population_dynamics(
#    change = growth(transition_matrix = egk_mat,
#                    transition_function = deterministic_transitions(transition_matrix)),
#    dispersal = NULL,
#    modification = NULL,
#    density_dependence = NULL)
#  

## ---- message = FALSE, progress = FALSE, fig.align = "center", out.width = "100%"----
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k)

egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        global_stochasticity = egk_mat_stoch),
                                        dispersal = NULL,
                                        modification = NULL,
                                        density_dependence = ceiling_density(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, fig.align = "center"------------------------------
plot(egk_results[1], type = "raster", stage = 2, timesteps = c(1, 10, 20), panels = c(3, 1))

## ---- message = FALSE, fig.align = "center"------------------------------
plot(egk_results[1], object = "carrying_capacity", type = "raster", timesteps = c(1, 10, 20), panels = c(3, 1))

## ---- message = FALSE, fig.align = "center", out.width = "100%"----------
egk_results[[2]][[5]][[1]][[1]]

## ---- message = FALSE, fig.align = "center", out.width = "100%"----------
egk_results[[2]][[3]][[2]]

## ---- message = FALSE, fig.align = "center", out.width = "100%"----------
egk_results[[2]][[3]][["carrying_capacity"]]

## ---- message = FALSE, fig.align = "center", out.width = "100%"----------
extract_spatial(egk_results, replicate = 2, timestep = 3, landscape_object = "carrying_capacity")

## ---- message = FALSE, fig.align = "center", out.width = "100%"----------
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(extract_spatial(egk_results), box = FALSE, axes = FALSE, col = viridis(100))

## ---- message = FALSE, progress = FALSE, fig.align = "center", out.width = "100%"----
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k,
                           source = egk_source,
                           sink = egk_sink)

egk_pop_dynamics <- population_dynamics(
  change = growth(transition_matrix = egk_mat,
                  global_stochasticity = egk_mat_stoch),
  modification = translocation(source_layer = "source",
                               sink_layer = "sink",
                               stages = 3,
                               effect_timesteps = c(1, 5, 10, 15)),
  density_dependence = NULL)

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- warning = FALSE, error = FALSE, message = FALSE, progress = FALSE, fig.align = "center", out.width = "100%"----
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k)

egk_pop_dynamics <- population_dynamics(
  change = growth(transition_matrix = egk_mat,
                  global_stochasticity = egk_mat_stoch),
  dispersal = kernel_dispersal(arrival_probability = "suitability",
                               max_distance = 5000,
                               dispersal_kernel = exponential_dispersal_kernel(distance_decay = 5000)),
  density_dependence = ceiling_density(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, fig.align = "center"------------------------------
plot(egk_results[1], type = "raster", stage = 3, timesteps = c(1, 10, 20), panels = c(3, 1))

## ---- message = FALSE, progress = FALSE, fig.align = "center", out.width = "100%"----
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k)

egk_pop_dynamics <- population_dynamics(
  change = growth(transition_matrix = egk_mat,
                  global_stochasticity = egk_mat_stoch),
  dispersal = kernel_dispersal(arrival_probability = "suitability",
                               max_distance = 5000,
                               dispersal_kernel = exponential_dispersal_kernel(distance_decay = 5000),
                               dispersal_proportion = carrying_capacity_dispersal()),
  modification = NULL,
  density_dependence = ceiling_density(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

