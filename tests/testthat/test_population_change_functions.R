context('population_change_functions-class')

test_that('population change functions class works', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k)
  
  two_sex_names <- c(paste0(colnames(egk_mat), "_F"), paste0(colnames(egk_mat), "_M"))
  
  egk_mat_2sex <- matrix(0, nrow = length(two_sex_names), ncol = length(two_sex_names))
  
  colnames(egk_mat_2sex) <- rownames(egk_mat_2sex) <- two_sex_names
  
  egk_mat_2sex[4:6, 4:6] <- egk_mat
  egk_mat_2sex[4, ] <- 0
  egk_mat_2sex[1:3, 1:3] <- egk_mat
  egk_mat_2sex[4, 1:3] <- egk_mat[1, ]
  
  egk_pop_2sex <- stack(egk_pop, egk_pop)
  names(egk_pop_2sex) <- two_sex_names
  
  landscape_2sex <- landscape(population = egk_pop_2sex,
                              suitability = egk_hab,
                              carrying_capacity = egk_k)

  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                 transition_order = "survival",
                                                 global_stochasticity = egk_mat_stoch * 0.5,
                                                 local_stochasticity = egk_mat_stoch * 0.5),
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL)
  
  pop_dyn_2sex <- population_dynamics(change = growth(transition_matrix = egk_mat_2sex,
                                                      two_sex = TRUE),
                                      dispersal = NULL,
                                      modification = NULL,
                                      density_dependence = NULL)
  
  pop_dyn_trans_fun <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                           transition_function = list(modified_transition(), competition_density())),
                                           dispersal = NULL,
                                           modification = NULL,
                                           density_dependence = NULL)
  
  pop_dyn_bad_trans_fun <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                           transition_function = list(1, competition_density())),
                                           dispersal = NULL,
                                           modification = NULL,
                                           density_dependence = NULL)
  
  pop_dyn_bad_mat_values <- population_dynamics(change = growth(transition_matrix = egk_mat * 2),
                                                dispersal = NULL,
                                                modification = NULL,
                                                density_dependence = NULL)
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim <- simulation(landscape = landscape_2sex,
                    population_dynamics = pop_dyn_2sex,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_trans_fun,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE,
                    demo_stochasticity = "none")
  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn_bad_mat_values,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
  
  expect_error(simulation(landscape = landscape,
                          population_dynamics = pop_dyn_bad_trans_fun,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE))
  
})
