context('visualisation-functions')

test_that('visualisation functions work', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k)
  
  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                 global_stochasticity = 0.05),
                                 dispersal = NULL,
                                 modification = NULL,
                                 density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                            population_dynamics = pop_dyn,
                            habitat_dynamics = NULL,
                            timesteps = 10,
                            replicates = 3,
                            verbose = FALSE)
  
plot_pop_trend(sim)

plot_k_trend(sim)

plot_pop_spatial(sim, timesteps = 1:9)

plot_k_spatial(sim, timesteps = 1:9)

plot_hab_spatial(sim, timesteps = 1:9)

compare_emp(sim, sim, all_points = TRUE)
  
})
