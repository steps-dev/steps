context('growth_transition_functions-class')

test_that('growth transition functions class works', {
  
  library(raster)
  library(future)
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k)
  
  pop_dyn <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                 transition_function = modified_transition(survival_layer = "suitability",
                                                                                           fecundity_layer = "suitability")),
                                 dispersal = cellular_automata_dispersal(),
                                 modification = NULL,
                                 density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  
  pop_dyn2 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                 transition_function = competition_density()),
                                 dispersal = cellular_automata_dispersal(),
                                 modification = NULL,
                                 density_dependence = NULL)
  
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn2,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  
})
