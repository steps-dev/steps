## ---- message = FALSE, echo = FALSE--------------------------------------
library(steps)
library(raster)

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

#for (i in seq_len(20)) {
#  writeRaster(egk_hab * 2, paste0("../working/rasters/egk_sf_", sprintf("%02i", i), "_0103.tif"), overwrite=TRUE)
#}

#for (i in seq_len(20)) {
#  writeRaster(egk_hab * 0.5, paste0("../working/rasters/egk_sf_", sprintf("%02i", i), "_0201.tif"), overwrite=TRUE)
#}

#for (i in seq_len(20)) {
#  writeRaster(egk_hab * 0.85, paste0("../working/rasters/egk_sf_", sprintf("%02i", i), "_0302.tif"), overwrite=TRUE)
#}

#for (i in seq_len(20)) {
#  writeRaster(egk_hab * 0.85, paste0("../working/rasters/egk_sf_", sprintf("%02i", i), "_0303.tif"), overwrite=TRUE)
#}


## ---- message = FALSE, fig.align="center"--------------------------------
egk_hab

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(egk_hab, box = FALSE, axes = FALSE)

## ---- message = FALSE, fig.align="center"--------------------------------
egk_k

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(egk_k, box = FALSE, axes = FALSE)

## ---- message = FALSE, fig.align="center"--------------------------------
egk_pop

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spplot(egk_pop)

## ---- message = FALSE----------------------------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = NULL, # could also specify suitability (egk_hab) here
                           carrying_capacity = NULL) # could also specify carrying capacity (egk_k) here

## ---- message = FALSE----------------------------------------------------
egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                        dispersal = NULL,
                                        modification = NULL,
                                        density_dependence = NULL)

## ---- message = FALSE,  results='hide', progress=FALSE-------------------
egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 1,
                          verbose = FALSE)

## ---- message = FALSE, fig.align="center"--------------------------------
plot(egk_results)

## ---- message = FALSE, fig.width=4, fig.align="center"-------------------
plot(egk_results, stage = 0)

## ---- message = FALSE, fig.width=4, fig.align="center"-------------------
plot(egk_results, stage = 2, newplot = TRUE)

## ---- message = FALSE, fig.align="center"--------------------------------
plot(egk_results, type = "raster", stage = 2)

## ---- message = FALSE, echo = FALSE--------------------------------------

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

## ---- message = FALSE, fig.width=4, fig.align="center"-------------------
plot(egk_results)

## ---- message = FALSE, fig.width=4, fig.align="center"-------------------
plot(egk_results[3], stage = 0)

## ---- message = FALSE, fig.align="center"--------------------------------
plot(egk_results[1], type = "raster", stage = 2)

## ---- message = FALSE----------------------------------------------------
egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        global_stochasticity = egk_mat_stoch),
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

## ---- message = FALSE----------------------------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = NULL)

egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
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

## ---- message = FALSE----------------------------------------------------
deterministic_transitions <- function(transition_matrix) {
  
  dim <- nrow(transition_matrix)
  
  function (landscape, timestep) {
    
    #### This assumes that the files are named with a particular convention
    #### and will not work for all cases. User must change code accordingly.
    
    #get metrics and constructor info
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
    current_timestep <- sprintf("%02i", timestep)
    n_cells <- length(which(!is.na(raster::getValues(landscape$population[[1]]))))
    
    #get relevant rasters - note, your working directory will be different
    files <- list.files("../working/rasters", pattern = paste0("_", current_timestep, "_")) 
    
    #initialise array
    transition_array <- array(0, dim = c(dim, dim, n_cells))
    
    #populate array:
    for (file in files) {
      r <- as.integer(substr(substr(file, nchar(file) - (8-1), nchar(file)), 1, 2))
      c <- as.integer(substr(substr(file, nchar(file) - (6-1), nchar(file)), 1, 2))
      #note, your working directory will be different
      transition_array[r, c, ] <- raster::raster(paste0("../working/rasters/", file))[cell_idx]
    }
    
    transition_array
  }
}
    
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = NULL)

egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        transition_function = deterministic_transitions(egk_mat)),
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

## ---- message = FALSE----------------------------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k)

egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        global_stochasticity = egk_mat_stoch),
                                        dispersal = NULL,
                                        modification = NULL,
                                        density_dependence = population_cap(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, fig.align="center"--------------------------------
plot(egk_results[1], object = "carrying_capacity")

## ---- message = FALSE, progress = FALSE----------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k)

egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        global_stochasticity = egk_mat_stoch),
                                        dispersal = kernel_dispersal(arrival_probability = "suitability"),
                                        modification = NULL,
                                        density_dependence = population_cap(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, fig.align="center"--------------------------------
plot(egk_results[1], type = "raster", stage = 3)

## ---- message = FALSE, progress = FALSE----------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k)

egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        global_stochasticity = egk_mat_stoch),
                                        dispersal = kernel_dispersal(arrival_probability = "suitability",
                                                                     dispersal_proportion = c(0, 0.5, 1)),
                                        modification = NULL,
                                        density_dependence = population_cap(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, progress = FALSE----------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k)

egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        global_stochasticity = egk_mat_stoch),
                                        dispersal = cellular_automata_dispersal(dispersal_distance = c(0, 5, 10),
                                                                                dispersal_proportion = c(0, 0.5, 1)),
                                        modification = NULL,
                                        density_dependence = population_cap(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, progress = FALSE----------------------------------
power_law_dispersal_kernel <- function () {
    fun <- function (r) 0.8*(1 + (r/0.9))^-2.1
  as.dispersal_function(fun)
}

egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k)

egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        global_stochasticity = egk_mat_stoch),
                                        dispersal = kernel_dispersal(arrival_probability = "suitability",
                                                                     dispersal_kernel = power_law_dispersal_kernel(),
                                                                     dispersal_proportion = c(0, 0.5, 1)),
                                        modification = NULL,
                                        density_dependence = population_cap(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE----------------------------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k,
                           fires = egk_dist)

egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                        global_stochasticity = egk_mat_stoch),
                                        dispersal = kernel_dispersal(arrival_probability = "suitability",
                                                                     dispersal_proportion = 0.5),
                                        modification = NULL,
                                        density_dependence = population_cap(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = list(disturbance(disturbance_layers = "fires",
                                                              effect_time = 3)),
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, fig.align="center"--------------------------------
plot(egk_results[1], type = "raster", stage = 3)

## ---- message = FALSE, fig.align="center"--------------------------------
plot(egk_results[1], object = "suitability")

