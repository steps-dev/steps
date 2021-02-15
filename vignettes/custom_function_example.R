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

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  
#  mortality <- function (mortality_layer, stages = NULL) {
#  
#    pop_dynamics <- function (landscape, timestep) {
#  
#        population_raster <- landscape$population
#        nstages <- raster::nlayers(population_raster)
#  
#        # get population as a matrix
#        idx <- which(!is.na(raster::getValues(population_raster[[1]])))
#        population_matrix <- raster::extract(population_raster, idx)
#  
#        mortality_prop <- raster::extract(landscape[[mortality_layer]][[timestep]], idx)
#  
#        if (is.null(stages)) stages <- seq_len(nstages)
#  
#        for (stage in stages) {
#  
#          population_matrix[ , stage] <- ceiling(population_matrix[ , stage] * mortality_prop)
#  
#        }
#  
#        # put back in the raster
#        population_raster[idx] <- population_matrix
#  
#        landscape$population <- population_raster
#  
#  
#      landscape
#  
#    }
#  
#    as.population_modification(pop_dynamics)
#  
#  }

## ---- message = FALSE, fig.align = "center"-----------------------------------
egk_pop

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spplot(egk_pop,
       col.regions = viridis(100),
       par.settings=list(strip.background = list(col = "white")))

## ---- message = FALSE---------------------------------------------------------
cellStats(egk_pop[[1]], sum)

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  percent_mortality <- function (percentage, stages = NULL) {
#  
#    pop_dynamics <- function (landscape, timestep) { ...

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  percent_mortality <- function (percentage, stages = NULL) {
#  
#    pop_dynamics <- function (landscape, timestep) {
#        population_raster <- landscape$population
#        nstages <- raster::nlayers(population_raster)
#  
#        # get population as a matrix
#        idx <- which(!is.na(raster::getValues(population_raster[[1]])))
#        population_matrix <- raster::extract(population_raster, idx) ...

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  percent_mortality <- function (percentage, stages = NULL) {
#  
#    pop_dynamics <- function (landscape, timestep) {
#        population_raster <- landscape$population
#        nstages <- raster::nlayers(population_raster)
#  
#        # get population as a matrix
#        idx <- which(!is.na(raster::getValues(population_raster[[1]])))
#        population_matrix <- raster::extract(population_raster, idx)
#  
#        if (is.null(stages)) stages <- seq_len(nstages)
#        for (stage in stages) {
#          initial_pop <- sum(population_matrix[ , stage])
#          changing_pop <- sum(population_matrix[ , stage])
#          while (changing_pop > initial_pop * percentage) {
#            non_zero <- which(population_matrix[idx , stage] > 0)
#            i <- sample(non_zero, 1)
#            population_matrix[i , stage] <- population_matrix[i , stage] - 1
#            changing_pop <- sum(population_matrix[ , stage])
#          }
#        } ...

## ---- message = FALSE---------------------------------------------------------
percent_mortality <- function (percentage, stages = NULL) {
  
  pop_dynamics <- function (landscape, timestep) {
    population_raster <- landscape$population
    nstages <- raster::nlayers(population_raster)
    
    # get population as a matrix
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population_matrix <- raster::extract(population_raster, idx)
    
    if (is.null(stages)) stages <- seq_len(nstages)
    for (stage in stages) {
      initial_pop <- sum(population_matrix[ , stage])
      changing_pop <- sum(population_matrix[ , stage])
      while (changing_pop > initial_pop * (1 - percentage)) {
        non_zero <- which(population_matrix[ , stage] > 0)
        i <- sample(non_zero, 1)
        population_matrix[i , stage] <- population_matrix[i , stage] - 1
        changing_pop <- sum(population_matrix[ , stage])
      }
    }
    
    population_raster[idx] <- population_matrix
    landscape$population <- population_raster
    landscape
  }
} 

## ---- message = FALSE---------------------------------------------------------
ls <- landscape(population = egk_pop,
                suitability = NULL,
                carrying_capacity = NULL)

pd <- population_dynamics(change = NULL,
                          dispersal = NULL,
                          modification = percent_mortality(percentage = 0.1,
                                                           stages = 1),
                          density_dependence = NULL)

results <- simulation(landscape = ls,
                      population_dynamics = pd,
                      habitat_dynamics = NULL,
                      demo_stochasticity = "none",
                      timesteps = 3,
                      replicates = 1,
                      verbose = FALSE)

## ---- message = FALSE---------------------------------------------------------
cbind(c("t0", "t1", "t2", "t3"),
      c(cellStats(egk_pop[[1]], sum),
        cellStats(results[[1]][[1]][["population"]][[1]], sum),
        cellStats(results[[1]][[2]][["population"]][[1]], sum),
        cellStats(results[[1]][[3]][["population"]][[1]], sum)))


## ---- message = FALSE, fig.align = "center"-----------------------------------
pop_stack <- stack(egk_pop[[1]],
                   results[[1]][[1]][["population"]][[1]],
                   results[[1]][[2]][["population"]][[1]],
                   results[[1]][[3]][["population"]][[1]])

names(pop_stack) <- c("t0", "t1", "t2", "t3")

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spplot(pop_stack,
       col.regions = viridis(100),
       par.settings=list(strip.background = list(col = "white")))

## ---- message = FALSE---------------------------------------------------------
percent_mortality_hab <- function (percentage, threshold, stages = NULL) {
  
  pop_dynamics <- function (landscape, timestep) {
    population_raster <- landscape$population
    nstages <- raster::nlayers(population_raster)
    
    # get population as a matrix
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population_matrix <- raster::extract(population_raster, idx)
    
    # get suitability cell indexes
    suitable <- raster::extract(landscape$suitability, idx)
    
    if (is.null(stages)) stages <- seq_len(nstages)
    for (stage in stages) {
      initial_pop <- sum(population_matrix[ , stage])
      changing_pop <- sum(population_matrix[ , stage])
      
      highly_suitable <- which(suitable >= threshold) # check for suitable cells
      
      while (changing_pop > initial_pop * (1 - percentage)) {
        non_zero <- which(population_matrix[ , stage] > 0)
        i <- sample(intersect(highly_suitable, non_zero), 1) # change the sampling pool
        population_matrix[i , stage] <- population_matrix[i , stage] - 1
        changing_pop <- sum(population_matrix[ , stage])
      }
    }
    
    population_raster[idx] <- population_matrix
    landscape$population <- population_raster
    landscape
  }
} 

## ---- message = FALSE---------------------------------------------------------
ls <- landscape(population = egk_pop,
                suitability = egk_hab,
                carrying_capacity = NULL)

pd <- population_dynamics(change = NULL,
                          dispersal = NULL,
                          modification = percent_mortality_hab(percentage = .10,
                                                               threshold = 0.7,
                                                               stages = 1),
                          density_dependence = NULL)

results <- simulation(landscape = ls,
                      population_dynamics = pd,
                      habitat_dynamics = NULL,
                      demo_stochasticity = "none",
                      timesteps = 3,
                      replicates = 1,
                      verbose = FALSE)

## ---- message = FALSE---------------------------------------------------------
cbind(c("t0", "t1", "t2", "t3"),
      c(cellStats(egk_pop[[1]], sum),
        cellStats(results[[1]][[1]][["population"]][[1]], sum),
        cellStats(results[[1]][[2]][["population"]][[1]], sum),
        cellStats(results[[1]][[3]][["population"]][[1]], sum)))


## ---- message = FALSE, fig.align = "center"-----------------------------------
pop_stack <- stack(egk_pop[[1]],
                   results[[1]][[1]][["population"]][[1]],
                   results[[1]][[2]][["population"]][[1]],
                   results[[1]][[3]][["population"]][[1]])

names(pop_stack) <- c("t0", "t1", "t2", "t3")

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spplot(pop_stack,
       col.regions = viridis(100),
       par.settings=list(strip.background = list(col = "white")))

## ---- message = FALSE, fig.align = "center"-----------------------------------
sampling_mask <- egk_hab # copy the habitat layer to assume its properties
sampling_mask[!is.na(egk_hab)] <- 0 # set all non-NA values to zero

# get index values of the top-right corner (15 x 15 cells)
idx_corner <- sort(unlist(lapply(0:15, function (x) seq(21, 541, by = 36) + x)))

# set non-NA raster values in top-right corner to one
sampling_mask[intersect(idx_corner, which(!is.na(egk_hab[])))] <- 1

plot(sampling_mask, axes = FALSE, box = FALSE)

ls <- landscape(population = egk_pop,
                suitability = NULL,
                carrying_capacity = NULL,
                "sampling_mask" = sampling_mask) # name the object to reference in our function

## ---- message = FALSE---------------------------------------------------------
percent_mortality_ras <- function (percentage, masking_raster, threshold, stages = NULL) {
  
  pop_dynamics <- function (landscape, timestep) {
    population_raster <- landscape$population
    nstages <- raster::nlayers(population_raster)
    
    # get population as a matrix
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population_matrix <- raster::extract(population_raster, idx)
    
    # Note, here we now use the name of the raster ('masking_raster' parameter)
    # to extract the object from the landscape object
    suitable <- raster::extract(landscape[[masking_raster]], idx)
    
    if (is.null(stages)) stages <- seq_len(nstages)
    for (stage in stages) {
      initial_pop <- sum(population_matrix[ , stage])
      changing_pop <- sum(population_matrix[ , stage])
      
      highly_suitable <- which(suitable >= threshold) # check for suitable cells
      
      while (changing_pop > initial_pop * (1 - percentage)) {
        non_zero <- which(population_matrix[ , stage] > 0)
        i <- sample(intersect(highly_suitable, non_zero), 1) # change the sampling pool
        population_matrix[i , stage] <- population_matrix[i , stage] - 1
        changing_pop <- sum(population_matrix[ , stage])
      }
    }
    
    population_raster[idx] <- population_matrix
    landscape$population <- population_raster
    landscape
  }
} 

## ---- message = FALSE---------------------------------------------------------
pd <- population_dynamics(change = NULL,
                          dispersal = NULL,
                          modification = percent_mortality_ras(percentage = .10,
                                                               masking_raster = "sampling_mask",
                                                               threshold = 0.7,
                                                               stages = 1),
                          density_dependence = NULL)

results <- simulation(landscape = ls,
                      population_dynamics = pd,
                      habitat_dynamics = NULL,
                      demo_stochasticity = "none",
                      timesteps = 3,
                      replicates = 1,
                      verbose = FALSE)

## ---- message = FALSE---------------------------------------------------------
cbind(c("t0", "t1", "t2", "t3"),
      c(cellStats(egk_pop[[1]], sum),
        cellStats(results[[1]][[1]][["population"]][[1]], sum),
        cellStats(results[[1]][[2]][["population"]][[1]], sum),
        cellStats(results[[1]][[3]][["population"]][[1]], sum)))


## ---- message = FALSE, fig.align = "center"-----------------------------------
pop_stack <- stack(egk_pop[[1]],
                   results[[1]][[1]][["population"]][[1]],
                   results[[1]][[2]][["population"]][[1]],
                   results[[1]][[3]][["population"]][[1]])

names(pop_stack) <- c("t0", "t1", "t2", "t3")

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spplot(pop_stack,
       col.regions = viridis(100),
       par.settings=list(strip.background = list(col = "white")))

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  modified_transition <- function(survival_layer = NULL,
#                                  fecundity_layer = NULL) {
#  
#    fun <- function (transition_array, landscape, timestep) {
#  
#      transition_matrix <- transition_array[, , 1]
#      idx <- which(transition_matrix != 0)
#      is_recruitment <- upper.tri(transition_matrix)[idx]
#  
#      array_length <- dim(transition_array)[3]
#  
#      cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
#  
#      if (is.null(survival_layer)) {
#        surv_mult <- rep(1, length(cell_idx))
#      } else {
#        if (raster::nlayers(landscape$suitability) > 1) {
#          surv_mult <- landscape[[survival_layer]][[timestep]][cell_idx]
#        } else {
#          surv_mult <- landscape[[survival_layer]][cell_idx]
#        }
#      }
#  
#      if (is.null(fecundity_layer)) {
#        fec_mult <- rep(1, length(cell_idx))
#      } else {
#        if (raster::nlayers(landscape$suitability) > 1) {
#          fec_mult <- landscape[[fecundity_layer]][[timestep]][cell_idx]
#        } else {
#          fec_mult <- landscape[[fecundity_layer]][cell_idx]
#        }
#      }
#  
#      for (i in seq_len(array_length)) {
#        new_surv <- transition_array[, , i][idx[!is_recruitment]] * surv_mult[i]
#        new_fec <- transition_array[, , i][idx[is_recruitment]] * fec_mult[i]
#  
#        transition_array[, , i][idx[!is_recruitment]] <- new_surv
#        transition_array[, , i][idx[is_recruitment]] <- new_fec
#      }
#  
#      transition_array
#  
#    }
#  
#    as.transition_function(fun)
#  
#  }

## ---- message = FALSE---------------------------------------------------------
ls <- landscape(population = egk_pop,
                suitability = NULL,
                carrying_capacity = NULL)

pd <- population_dynamics(change = growth(transition_matrix = egk_mat),
                          dispersal = NULL,
                          modification = NULL,
                          density_dependence = NULL)

results <- simulation(landscape = ls,
                      population_dynamics = pd,
                      habitat_dynamics = NULL,
                      timesteps = 5,
                      replicates = 1,
                      verbose = FALSE)

## ---- message = FALSE, fig.align = "center"-----------------------------------
plot(results, stages = 1)

## ---- message = FALSE, fig.align = "center"-----------------------------------
plot(results,
     stages = 1,
     type = "raster",
     timesteps = 1:5,
     panels = c(5, 1))

## ---- message = FALSE, echo = FALSE-------------------------------------------
print(egk_mat)

paste("Rmax =", abs(eigen(egk_mat)$values[1]))

## ---- message = FALSE, fig.align = "center", warning = FALSE------------------
adult_fec <- egk_hab # copy the habitat layer to assume its properties
adult_fec[] <- NA # set all values to NA

ncells <- ncell(adult_fec)

# assign random numbers between 0.1 and 0.4 to 80% of the cells
adult_fec[sample(seq_len(ncells), 0.8 * ncells, replace = FALSE)] <- runif(ncells, 0.1, 0.5)

plot(adult_fec, axes = FALSE, box = FALSE)

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  define_fecundity <- function(fecundity_layer) {
#  
#    fun <- function (transition_array, landscape, timestep) { ...

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  define_fecundity <- function(fecundity_layer) {
#  
#    fun <- function (transition_array, landscape, timestep) {
#      transition_matrix <- transition_array[, , 1]
#      idx <- which(transition_matrix != 0)
#      is_recruitment <- upper.tri(transition_matrix)[idx]
#      cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])
#  
#      # this length should be pre-defined by non-NA cells in the landscape
#      array_length <- dim(transition_array)[3]
#      ...

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  define_fecundity <- function(fecundity_layer) {
#  
#    fun <- function (transition_array, landscape, timestep) {
#      transition_matrix <- transition_array[, , 1]
#      idx <- which(transition_matrix != 0)
#      is_recruitment <- upper.tri(transition_matrix)[idx]
#      cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])
#  
#      # this length should be pre-defined by non-NA cells in the landscape
#      array_length <- dim(transition_array)[3]
#  
#      fec_values <- landscape[[fecundity_layer]][cell_idx]
#  
#      for (i in seq_len(array_length)) {
#        new_fec <- fec_values[i]
#        if (!is.na(new_fec)) {
#          transition_array[, , i][tail(idx[is_recruitment], 1)] <- new_fec
#        }
#      } ...

## ---- message = FALSE---------------------------------------------------------
define_fecundity <- function(fecundity_layer) {
  
  fun <- function (transition_array, landscape, timestep) {
    transition_matrix <- transition_array[, , 1]
    idx <- which(transition_matrix != 0)
    is_recruitment <- upper.tri(transition_matrix)[idx]
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
    
    # this length should be pre-defined by non-NA cells in the landscape
    array_length <- dim(transition_array)[3]
          
    fec_values <- landscape[[fecundity_layer]][cell_idx]
    
    for (i in seq_len(array_length)) {
      new_fec <- fec_values[i]
      if (!is.na(new_fec)) {
        transition_array[, , i][tail(idx[is_recruitment], 1)] <- new_fec
      }
    }

    transition_array
  }
}

## ---- message = FALSE---------------------------------------------------------
ls <- landscape(population = egk_pop,
                suitability = NULL,
                carrying_capacity = NULL,
                "adult_fec" = adult_fec)

pd <- population_dynamics(
  change = growth(transition_matrix = egk_mat,
                  transition_function = define_fecundity(fecundity_layer = "adult_fec")),
  dispersal = NULL,
  modification = NULL,
  density_dependence = NULL)

results <- simulation(landscape = ls,
                      population_dynamics = pd,
                      habitat_dynamics = NULL,
                      timesteps = 5,
                      replicates = 1,
                      verbose = FALSE)

## ---- message = FALSE, fig.align = "center"-----------------------------------
plot(results, stages = 1)

## ---- message = FALSE, fig.align = "center"-----------------------------------
plot(results,
     stages = 1,
     type = "raster",
     timesteps = 1:5,
     panels = c(5, 1))

## ---- message = FALSE, eval = FALSE-------------------------------------------
#  exponential_dispersal_kernel <- function (distance_decay = 0.5, normalize = FALSE) {
#    if (normalize) {
#      fun <- function (r) (1 / (2 * pi * distance_decay ^ 2)) * exp(-r / distance_decay)
#    } else {
#      fun <- function (r) exp(-r / distance_decay)
#    }
#    as.dispersal_function(fun)
#  }

## ---- message = FALSE---------------------------------------------------------
r <- seq(1, 1000, 10)
distance_decay = 100
plot(r, exp(-r / distance_decay), type = 'l')

## ---- message = FALSE---------------------------------------------------------
logistic_dispersal_kernel <- function (dist_shape, dist_slope) {
fun <- function (dist) 1 / (1 + exp((dist - dist_shape) / dist_slope))
fun
}

r <- seq(1, 1000, 10)
dist_shape = 500
dist_slope = 20
plot(r, 1 / (1 + exp((r - dist_shape) / dist_slope)), type = 'l')

## ---- message = FALSE---------------------------------------------------------
egk_pop_thin <- egk_pop
egk_pop_thin[sample(1:ncell(egk_pop), 0.99 * ncell(egk_pop))] <- 0

ls <- landscape(population = egk_pop_thin,
                suitability = NULL,
                carrying_capacity = NULL)

pd_exp_disp <- population_dynamics(
  change = NULL,
  dispersal = kernel_dispersal(dispersal_kernel = exponential_dispersal_kernel(distance_decay = 100),
                               max_distance = 1000,
                               arrival_probability = "none"),
  modification = NULL,
  density_dependence = NULL)

results_01 <- simulation(landscape = ls,
                      population_dynamics = pd_exp_disp,
                      habitat_dynamics = NULL,
                      timesteps = 10,
                      replicates = 1,
                      verbose = FALSE)

plot(results_01, type = "raster", stages = 3, timesteps = c(1,5,10), panels = c(3, 1))

pd_log_disp <- population_dynamics(
  change = NULL,
  dispersal = kernel_dispersal(dispersal_kernel = logistic_dispersal_kernel(dist_shape = 500, dist_slope = 20),
                               max_distance = 1000,
                               arrival_probability = "none"),
  modification = NULL,
  density_dependence = NULL)

results_02 <- simulation(landscape = ls,
                      population_dynamics = pd_log_disp,
                      habitat_dynamics = NULL,
                      timesteps = 10,
                      replicates = 1,
                      verbose = FALSE)

plot(results_02, type = "raster", stages = 3, timesteps = c(1,5,10), panels = c(3, 1))

## ---- message = FALSE---------------------------------------------------------
plan(multisession)
results_02 <- simulation(landscape = ls,
                      population_dynamics = pd_log_disp,
                      habitat_dynamics = NULL,
                      timesteps = 10,
                      replicates = 3,
                      verbose = FALSE,
                      future.globals = list("logistic_dispersal_kernel" = logistic_dispersal_kernel))

