context('population_change_functions-class')

test_that('population change functions class works', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k)
  
  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                 transition_order = "survival",
                                                 global_stochasticity = egk_mat_stoch,
                                                 local_stochasticity = egk_mat_stoch),
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  
  two_sex_names <- c(paste0(colnames(egk_mat), "_F"), paste0(colnames(egk_mat), "_M"))
  
  egk_mat2 <- matrix(0, nrow = length(two_sex_names), ncol = length(two_sex_names))
  
  colnames(egk_mat2) <- rownames(egk_mat2) <- two_sex_names
  
  egk_mat2[4:6, 4:6] <- egk_mat
  egk_mat2[4, ] <- 0
  egk_mat2[1:3, 1:3] <- egk_mat
  egk_mat2[4, 1:3] <- egk_mat[1, ]
  
  egk_pop2 <- stack(egk_pop, egk_pop)
  names(egk_pop2) <- two_sex_names
  
  landscape2 <- landscape(population = egk_pop2,
                          suitability = egk_hab,
                          carrying_capacity = egk_k)
  
  
  pop_dyn2 <- population_dynamics(change = growth(transition_matrix = egk_mat2,
                                                 two_sex = TRUE),
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape2,
                    population_dynamics = pop_dyn2,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  
})
