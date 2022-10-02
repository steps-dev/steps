context('visualisation-functions')

test_that('visualisation functions work', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k)
  
  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                 global_stochasticity = 0.01),
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL)
  
  pop_dyn_crash <- population_dynamics(change = growth(transition_matrix = egk_mat * 0.5,
                                                       global_stochasticity = 0.005),
                                       dispersal = NULL,
                                       modification = NULL,
                                       density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim_crash <- simulation(landscape = landscape,
                          population_dynamics = pop_dyn_crash,
                          habitat_dynamics = NULL,
                          timesteps = 10,
                          replicates = 3,
                          verbose = FALSE)
  
  plot_pop_trend(sim)
  
  plot_pop_trend(sim_crash, emp = TRUE)
  
  plot_k_trend(sim)
  
  out <- plot_k_trend(sim, summary_stat = "sum", return_data = TRUE)
  
  plot_pop_spatial(sim)
  
  plot_pop_spatial(sim, stage = 2, timesteps = 1:9)
  
  plot_k_spatial(sim)
  
  plot_k_spatial(sim, timesteps = 1:9)
  
  plot_hab_spatial(sim)
  
  plot_hab_spatial(sim, timesteps = 1:9)
  
  compare_emp(sim, sim_crash, all_points = TRUE)
  
})
