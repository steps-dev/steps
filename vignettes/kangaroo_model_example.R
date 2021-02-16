## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  dpi = 300,
  fig.width = 7,
  out.width = "100%",
  cache = TRUE
)

## ---- message = FALSE---------------------------------------------------------
library(steps)
library(raster)
library(viridisLite)
library(future)

## ---- message = FALSE---------------------------------------------------------

egk_mat <- matrix(c(0.00,0.55,0.85,
                    0.50,0.20,0.00,
                    0.00,0.55,0.85),
                  nrow = 3,
                  ncol = 3,
                  byrow = TRUE)
colnames(egk_mat) <- rownames(egk_mat) <- c('juvenile','subadult','adult')

egk_mat_stoch <- matrix(c(0.00,0.05,0.05,
                          0.05,0.05,0.00,
                          0.00,0.01,0.01),
                        nrow = 3,
                        ncol = 3,
                        byrow = TRUE)
colnames(egk_mat_stoch) <- rownames(egk_mat_stoch) <- c('juvenile','subadult','adult')

## ---- message = FALSE, fig.align = "center"-----------------------------------
egk_hab

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(egk_hab, box = FALSE, axes = FALSE, col = viridis(100), main = "Habitat Suitability")

## ---- message = FALSE, fig.align = "center"-----------------------------------
egk_k

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(egk_k, box = FALSE, axes = FALSE, col = viridis(100))

## ---- message = FALSE, fig.align = "center"-----------------------------------
egk_pop

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spplot(egk_pop, col.regions = viridis(100))

## ---- message = FALSE---------------------------------------------------------
library(foreach)
stable_states <- abs( eigen(egk_mat)$vectors[,1] / sum(eigen(egk_mat)$vectors[,1]))
popN <- stack(replicate(ncol(egk_mat), egk_k)) * stable_states
idx <- which(!is.na(getValues(popN[[1]])))
pop <- stack(
  foreach(i = 1:nlayers(popN)) %do% {
    max_pop <- ceiling(cellStats(popN[[i]], max, na.rm = T))
    pop_values <- popN[[i]][idx]
    popN[[i]][idx] <- rbinom(prob = (pop_values/max_pop),
                             size = max_pop,
                             n = length(pop_values))
    popN[[i]]
  })
names(pop) <- colnames(egk_mat)
pop

## ---- message = FALSE---------------------------------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = NULL,
                           carrying_capacity = NULL)

## ---- message = FALSE---------------------------------------------------------
egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat),
                                        dispersal = NULL,
                                        modification = NULL,
                                        density_dependence = NULL)

## ---- message = FALSE,  results = 'hide', progress = FALSE--------------------
egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 1,
                          verbose = FALSE)

## ---- message = FALSE, fig.align = "center"-----------------------------------
plot(egk_results)

## ---- message = FALSE, fig.align = "center", fig.width = 3, fig.height = 4, out.width = "33%"----
plot(egk_results, stage = 0)

## ---- message = FALSE, fig.align = "center", fig.width = 3, fig.height = 4, out.width = "33%"----
plot(egk_results, stage = 2)

## ---- message = FALSE, fig.align = "center", fig.show='hold'------------------
plot(egk_results, type = "raster", stage = 2, timesteps = c(1, 10, 20), panels = c(3, 1))

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  plot(egk_results, type = "raster", stage = 2, timesteps = c(1, 10, 20), animate = TRUE)

## ---- message = FALSE, progress = FALSE, eval = FALSE-------------------------
#  plan(multisession, workers = 3) # This is how we specify to simulate
#  # replicates on separate processors in parallel
#  
#  egk_results <- simulation(landscape = egk_landscape,
#                            population_dynamics = egk_pop_dynamics,
#                            habitat_dynamics = NULL,
#                            timesteps = 20,
#                            replicates = 3,
#                            verbose = FALSE)
#  

## ---- message = FALSE, fig.align = "center"-----------------------------------
plot(egk_results[1], type = "raster", stage = 0, timesteps = c(1, 10, 20), panels = c(3, 1))

## ---- message = FALSE, progress = FALSE, fig.align = "center"-----------------
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

## ---- message = FALSE, progress = FALSE, fig.align = "center"-----------------
plan(sequential)

egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = NULL)

egk_pop_dynamics <- population_dynamics(
  change = growth(transition_matrix = egk_mat,
                  global_stochasticity = egk_mat_stoch,
                  transition_function = modified_transition(survival_layer = "suitability",
                                                            fecundity_layer = NULL)),
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

## ---- message = FALSE, progress = FALSE, fig.align = "center"-----------------
deterministic_transitions <- function(spatial_object) {
  
  fun <- function (transition_array, landscape, timestep) {
    
    #### This assumes that the spatial layers in the raster stack are named with
    #### a particular convention and will not work for all cases.
    
    # get metrics and constructor info
    current_timestep <- sprintf("%02i", timestep)
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
    n_cells <- length(which(!is.na(raster::getValues(landscape$population[[1]]))))
    
    # get names of rasters and subset by timestep
    names_all <- names(landscape[[spatial_object]])
    names_timestep <- grep(paste0("^.*", current_timestep, "\\_"), names_all, value = TRUE)
    
    # populate array:
    for (name in names_timestep) {
      r <- as.integer(substr(substr(name, nchar(name) - 3, nchar(name)), 1, 2))
      c <- as.integer(substr(substr(name, nchar(name) - 1, nchar(name)), 1, 2))
      transition_array[r, c, ] <- landscape[[spatial_object]][[name]][cell_idx]
    }
    
    #### Return array with required dimensions
    transition_array
  }

}

egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = NULL,
                           "sf_mods" = egk_sf)

egk_pop_dynamics <- population_dynamics(
  change = growth(transition_matrix = egk_mat,
                  transition_function = deterministic_transitions(spatial_object = "sf_mods")),
  dispersal = NULL,
  modification = NULL,
  density_dependence = NULL)

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, progress = FALSE, fig.align = "center", eval = FALSE----
#  egk_landscape <- landscape(population = egk_pop,
#                             suitability = egk_hab,
#                             carrying_capacity = egk_k)
#  
#  egk_pop_dynamics <- population_dynamics(change = growth(transition_matrix = egk_mat,
#                                                          global_stochasticity = egk_mat_stoch),
#                                          dispersal = NULL,
#                                          modification = NULL,
#                                          density_dependence = ceiling_density(stages = 3))
#  
#  egk_results <- simulation(landscape = egk_landscape,
#                            population_dynamics = egk_pop_dynamics,
#                            habitat_dynamics = NULL,
#                            timesteps = 20,
#                            replicates = 3,
#                            verbose = FALSE)

## ---- message = FALSE, progress = FALSE, fig.align = "center", eval = FALSE----
#  egk_landscape <- landscape(population = egk_pop,
#                             suitability = NULL,
#                             carrying_capacity = NULL,
#                             "origins" = egk_origins,
#                             "destinations" = egk_destinations)
#  
#  egk_pop_dynamics <- population_dynamics(
#    change = growth(transition_matrix = egk_mat,
#                    global_stochasticity = egk_mat_stoch),
#    dispersal = NULL,
#    modification = translocation(origins_layer = "origins",
#                                 destinations_layer = "destinations",
#                                 stages = 3,
#                                 effect_timesteps = c(1, 5, 10, 15)),
#    density_dependence = NULL)
#  
#  egk_results <- simulation(landscape = egk_landscape,
#                            population_dynamics = egk_pop_dynamics,
#                            habitat_dynamics = NULL,
#                            timesteps = 20,
#                            replicates = 3,
#                            verbose = FALSE)

## ---- warning = FALSE, error = FALSE, message = FALSE, progress = FALSE, fig.align = "center"----
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k)

egk_pop_dynamics <- population_dynamics(
  change = growth(transition_matrix = egk_mat,
                  global_stochasticity = egk_mat_stoch),
  dispersal = kernel_dispersal(
    arrival_probability = "suitability",
    max_distance = 5000,
    dispersal_kernel = exponential_dispersal_kernel(distance_decay = 5000)
  ),
  density_dependence = ceiling_density(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

## ---- message = FALSE, fig.align = "center"-----------------------------------
plot(egk_results[1], type = "raster", stage = 3, timesteps = c(1, 10, 20), panels = c(3, 1))

## ---- message = FALSE, progress = FALSE, fig.align = "center", eval = FALSE----
#  egk_landscape <- landscape(population = egk_pop,
#                             suitability = egk_hab,
#                             carrying_capacity = egk_k)
#  
#  egk_pop_dynamics <- population_dynamics(
#    change = growth(transition_matrix = egk_mat,
#                    global_stochasticity = egk_mat_stoch),
#    dispersal = kernel_dispersal(
#      arrival_probability = "suitability",
#      max_distance = 5000,
#      dispersal_kernel = exponential_dispersal_kernel(distance_decay = 5000),
#      dispersal_proportion = density_dependence_dispersing()
#    ),
#    modification = NULL,
#    density_dependence = ceiling_density(stages = 3))
#  
#  egk_results <- simulation(landscape = egk_landscape,
#                            population_dynamics = egk_pop_dynamics,
#                            habitat_dynamics = NULL,
#                            timesteps = 20,
#                            replicates = 3,
#                            verbose = FALSE)

## ---- message = FALSE, progress = FALSE, fig.align = "center", eval = FALSE----
#  egk_landscape <- landscape(population = egk_pop,
#                             suitability = egk_hab,
#                             carrying_capacity = egk_k)
#  
#  egk_pop_dynamics <- population_dynamics(
#    change = growth(transition_matrix = egk_mat,
#                    global_stochasticity = egk_mat_stoch),
#    dispersal = kernel_dispersal(arrival_probability = "suitability",
#                                 max_distance = 5000,
#                                 dispersal_kernel = function (r) exp(-r / 2000),
#                                 dispersal_proportion = density_dependence_dispersing()),
#    modification = NULL,
#    density_dependence = ceiling_density(stages = 3))
#  
#  egk_results <- simulation(landscape = egk_landscape,
#                            population_dynamics = egk_pop_dynamics,
#                            habitat_dynamics = NULL,
#                            timesteps = 20,
#                            replicates = 3,
#                            verbose = FALSE)

## ---- message = FALSE, progress = FALSE, fig.align = "center", eval = FALSE----
#  egk_landscape <- landscape(population = egk_pop,
#                             suitability = egk_hab,
#                             carrying_capacity = egk_k)
#  
#  egk_pop_dynamics <- population_dynamics(
#    change = growth(transition_matrix = egk_mat,
#                    global_stochasticity = egk_mat_stoch),
#    dispersal = cellular_automata_dispersal(
#      max_cells = c(0, 10, 20)
#    ),
#    modification = NULL,
#    density_dependence = ceiling_density(stages = 3))
#  
#  egk_results <- simulation(landscape = egk_landscape,
#                            population_dynamics = egk_pop_dynamics,
#                            habitat_dynamics = NULL,
#                            timesteps = 20,
#                            replicates = 3,
#                            verbose = FALSE)

## ---- message = FALSE, results='hide'-----------------------------------------
egk_landscape <- landscape(population = egk_pop,
                           suitability = egk_hab,
                           carrying_capacity = egk_k,
                           "fires" = egk_fire)

egk_pop_dynamics <- population_dynamics(
  change = growth(transition_matrix = egk_mat,
                  global_stochasticity = egk_mat_stoch),
  dispersal = kernel_dispersal(
    arrival_probability = "suitability",
    max_distance = 5000,
    dispersal_kernel = exponential_dispersal_kernel(distance_decay = 5000)
  ),
  modification = NULL,
  density_dependence = ceiling_density(stages = 3))

egk_results <- simulation(landscape = egk_landscape,
                          population_dynamics = egk_pop_dynamics,
                          habitat_dynamics = list(
                            fire_effects(fire_layers = "fires",
                                         effect_time = 5,
                                         regeneration_function = function (time) {-time})
                          ),
                          timesteps = 20,
                          replicates = 3,
                          verbose = FALSE)

plot(egk_results)

## ---- message = FALSE, fig.align = "center"-----------------------------------
plot(egk_results[1], type = "raster", stage = 3, timesteps = c(1, 10, 20), panels = c(3, 1))

## ---- message = FALSE, fig.align = "center"-----------------------------------
plot(egk_results[1], object = "suitability", timesteps = c(1, 10, 20), panels = c(3, 1))

## ---- message = FALSE, results='hide', eval = FALSE---------------------------
#  carrying_cap_fun <- function (landscape, timestep) {
#  
#    fun <- function(suitability) {
#      75 - round(75 * dlogis(suitability, scale = 0.25))
#    }
#  
#    suit <- landscape$suitability
#    if (raster::nlayers(suit) > 1) {
#      suit <- suit[[timestep]]
#    }
#  
#    raster::calc(suit, fun)
#  
#  }
#  
#  suit_seq <- seq(0, 1, 0.1)
#  plot(suit_seq, 75 - round(75 * dlogis(suit_seq, scale = 0.25)), type = 'l')
#  
#  
#  egk_landscape <- landscape(population = egk_pop,
#                             suitability = egk_hab,
#                             carrying_capacity = carrying_cap_fun,
#                             fires = egk_fire)
#  
#  egk_pop_dynamics <- population_dynamics(
#    change = growth(transition_matrix = egk_mat,
#                    global_stochasticity = egk_mat_stoch),
#    dispersal = kernel_dispersal(
#      arrival_probability = "suitability",
#      max_distance = 5000,
#      dispersal_kernel = exponential_dispersal_kernel(distance_decay = 5000)
#    ),
#    modification = NULL,
#    density_dependence = ceiling_density(stages = 3))
#  
#  egk_results <- simulation(landscape = egk_landscape,
#                            population_dynamics = egk_pop_dynamics,
#                            habitat_dynamics = list(fire_effects(fire_layers = "fires",
#                                                                 effect_time = 5,
#                                                                 regeneration_function =
#                                                                   function (time) {-time})),
#                            timesteps = 20,
#                            replicates = 3,
#                            verbose = FALSE)

## ---- message = FALSE, results = 'hide', eval = FALSE-------------------------
#  egk_landscape <- landscape(population = egk_pop,
#                             suitability = egk_hab,
#                             carrying_capacity = carrying_cap_fun,
#                             "fires" = egk_fire,
#                             "roads" = egk_road)
#  
#  egk_pop_dynamics <- population_dynamics(
#    change = growth(transition_matrix = egk_mat,
#                    global_stochasticity = egk_mat_stoch),
#    dispersal = kernel_dispersal(
#      arrival_probability = "suitability",
#      max_distance = 5000,
#      dispersal_kernel = exponential_dispersal_kernel(distance_decay = 5000)
#    ),
#    modification = NULL,
#    density_dependence = ceiling_density(stages = 3))
#  
#  egk_results <- simulation(landscape = egk_landscape,
#                            population_dynamics = egk_pop_dynamics,
#                            habitat_dynamics = list(
#                              fire_effects(fire_layers = "fires",
#                                           effect_time = 5,
#                                           regeneration_function = function (time) {-time}),
#                              disturbance(disturbance_layers = "roads",
#                                          effect_time = 1)
#                            ),
#                            timesteps = 20,
#                            replicates = 3,
#                            verbose = FALSE)

## ---- message = FALSE, fig.align = "center"-----------------------------------
egk_results[[2]][[5]][[1]][[1]]

## ---- message = FALSE, fig.align = "center"-----------------------------------
egk_results[[2]][[3]][[2]]

## ---- message = FALSE, fig.align = "center"-----------------------------------
egk_results[[2]][[3]][["suitability"]]

## ---- message = FALSE, fig.align = "center"-----------------------------------
extract_spatial(egk_results, replicate = 2, timestep = 3, landscape_object = "suitability")

## ---- message = FALSE, fig.align = "center"-----------------------------------
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(extract_spatial(egk_results), box = FALSE, axes = FALSE, col = viridis(100))

