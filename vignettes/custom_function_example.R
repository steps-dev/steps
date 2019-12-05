## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  dpi = 300,
  fig.width = 7,
  out.width = "100%",
  cache = TRUE
)

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

