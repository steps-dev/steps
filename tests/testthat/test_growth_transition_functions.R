context('growth_transition_functions-class')

test_that('growth transition functions class works', {
  
  library(raster)
  library(future)
  
  egk_mat_mask <- egk_mat > 0
  
  egk_hab_stack <- stack(replicate(10, egk_hab))
  
  landscape <- landscape(population = egk_pop,
                         suitability = egk_hab,
                         carrying_capacity = egk_k)
  
  landscape2 <- landscape(population = egk_pop,
                          suitability = egk_hab_stack,
                          carrying_capacity = egk_k)
  
  pop_dyn_mod <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                     transition_function = modified_transition(survival_layer = "suitability",
                                                                                               fecundity_layer = "suitability")),
                                     dispersal = cellular_automata_dispersal(),
                                     modification = NULL,
                                     density_dependence = NULL)

  pop_dyn_comp <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                      transition_function = competition_density()),
                                      dispersal = cellular_automata_dispersal(),
                                      modification = NULL,
                                      density_dependence = NULL)
  
  pop_dyn_comp2 <- population_dynamics(change = growth(transition_matrix = egk_mat,
                                                       transition_function = competition_density(mask = egk_mat_mask, stages = 2)),
                                       dispersal = cellular_automata_dispersal(),
                                       modification = NULL,
                                       density_dependence = NULL)
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_mod,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim <- simulation(landscape = landscape2,
                    population_dynamics = pop_dyn_mod,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_comp,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  sim <- simulation(landscape = landscape,
                    population_dynamics = pop_dyn_comp2,
                    habitat_dynamics = NULL,
                    timesteps = 10,
                    replicates = 3,
                    verbose = FALSE)
  
  
})
