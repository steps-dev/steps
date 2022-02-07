# skeleton architecture for dhmpr 

# ~~~~~~~~~~~~
# functions to make fake data for examples (not needed in package)

# fake a demographic transition matrix of the right dimensions
fake_transition_matrix <- function (n_stages) {
  
  survival <- runif(n_stages, 0.7, 0.9)
  growth <- runif(n_stages - 1, 0.5, 0.7)
  recruitment <- rlnorm(1)
  
  # base matrix 
  transition_matrix <- diag(n_stages)
  
  # add growth
  growth_idx <- which(row(transition_matrix) == col(transition_matrix) + 1)
  transition_matrix[growth_idx] <- growth
  
  # columns sum to 1
  sums <- colSums(transition_matrix)
  transition_matrix <- sweep(transition_matrix, 2, sums, FUN = "/")
  
  # apply survival
  transition_matrix <- sweep(transition_matrix, 2, survival, FUN = "*")
  
  # add recruitment
  transition_matrix[1, n_stages] <- recruitment
  
  transition_matrix
}

# make a multistage population matrix of the right dimensions (assuming each
# stage contributes equally to the carrying capacity)
fake_population <- function (carrying_capacity, n_stages) {
  
  # population of all stages
  total_population <- rbinom(n_patches, carrying_capacity, 0.5)
  
  # stable stage distribution
  ssd <- rev(cumsum(runif(n_stages)))
  ssd <- ssd / sum(ssd)
  
  # stage-specific populations
  population <- vapply(total_population,
                       function (N) {rmultinom(1, N, ssd)},
                       1:3)
  
  population <- t(population)
  
  population
  
}

# ~~~~~~~~~~~~
# internal functions - not to be exported

#  simple dispersal probability matrix from patch locations
dispersal_matrix <- function (locations, distance_decay = 0.5) {
  D <- as.matrix(dist(locations))
  dispersal_matrix <- exp(D / distance_decay)
  sums <- colSums(dispersal_matrix)
  dispersal_matrix <- sweep(dispersal_matrix, 2, sums, "/")
  dispersal_matrix
}

# for each timestep, do all the functions, in the order provided
iterate_system <- function (state, timesteps, setup) {
  
  output_states <- list()
  
  for (timestep in timesteps) {
    for (fun in setup) {
      state <- fun(state, timestep)
    }
   output_states[[timestep]] <- state
  }
  
  output_states
  
}

cap_population <- function (new_population, state) {
  
  carrying_capacity <- state$habitat$carrying_capacity
  
  if (is.null(carrying_capacity)) {
    stop ("carrying capacity must be specified",
          call. = FALSE)
  }
  
  # get degree of overpopulation, and shrink accordingly
  overpopulation <- carrying_capacity / rowSums(new_population)
  overpopulation <- pmin(overpopulation, 1)
  new_population <- sweep(new_population, 1, overpopulation, "*")
  
  new_population
}

check_dynamics_order <- function (order) {
  sorted_order <- sort(order)
  expected <- c("demographic_dynamics",
                "dispersal",
                "habitat_dynamics",
                "transition")
  if (!identical(sorted_order, expected)) {
    msg <- paste0("order must be a length-4 character vector giving the order ",
                  "in which to run the dynamic functions. It must contain each ",
                  "of the following strings once and only once:\n",
                  "'", paste(expected, collapse = "', '"), "'")
    stop (msg, call. = FALSE)
  }
}

set_class <- function (x, class) {
  class(x) <- c(class, class(x))
  x
}

as_dynamics <- function (object, type) {
  object <- set_class(object, "dynamics")
  object <- set_class(object, type)
  object
}

print.dynamics <- function (x, ...) {
  cat("this is a dynamics object")
}

print.dispersal_dynamics <- function (x, ...) {
  cat("this is a dispersal dynamics object")
}

# ~~~~~~~~~~~~
# exported package functions

# default *transition functions* defining how the system changes over time

# simple demographic transition (projection of fractional populations)
simple_transition <- function (state, timestep) {
  
  population <- state$habitat$population
  transition_matrix <- state$demography$transition_matrix
  
  new_population <- population %*% transition_matrix
  
  state$habitat$population <- new_population
  state
  
}

simple_transition <- as_transition(simple_transition, "demographic")

# simple dispersal (fractional)
simple_dispersal <- function (dispersal_decay) {

  dispersal <- function (state, timestep) {
    
    population <- state$habitat$population
    dispersal <- dispersal_matrix(state$habitat$locations)
    
    new_population <- dispersal %*% population
    
    state$habitat$population <- new_population
    state
    
  }
  
  as_transition(dispersal, "dispersal")
  
}
  

# no changes to the habitat
no_habitat_dynamics <- function (state, timestep) {
  state
}
no_habitat_dynamics <- as_transition(no_habitat_dynamics, "habitat_dynamics")

# no changes to the demographics
no_demographic_dynamics <- function (state, timestep) {
  state
}
no_demographic_dynamics <- as_transition(no_demographic_dynamics, "demographic_dynamics")

# some other options for transition functions
# demographic transition including carrying capacity
simple_transition_capped <- function (state, timestep) {
  
  population <- state$habitat$population
  transition_matrix <- state$demography$transition_matrix

  new_population <- population %*% transition_matrix
  
  new_population <- cap_population(new_population, state)
  
  state$habitat$population <- new_population
  state
  
}

# simple dispersal (fractional)
simple_dispersal_capped <- function (state, timestep) {
  
  population <- state$habitat$population
  dispersal <- dispersal_matrix(state$habitat$locations)
  
  new_population <- dispersal %*% population
  
  new_population <- cap_population(new_population, state)
  
  state$habitat$population <- new_population
  state
  
}

# stachasticly alter the demographic transition matrix, about a global value
environmental_stochasticity <- function (global_transition_matrix, stochasticity) {
  dim <- nrow(global_transition_matrix)
  idx <- which(global_transition_matrix != 0)
  recruitment_mask <- idx == ((dim ^ 2) - dim + 1)
  lower <- 0
  upper <- ifelse(recruitment_mask, 1, Inf)
  vals <- global_transition_matrix[idx]
  
  demographic_dynamics <- function (state, timestep) {
    
    transition_matrix <- global_transition_matrix
    
    
    transition_matrix[idx] <- extraDistr::rtnorm(length(idx),
                                                 vals,
                                                 stochasticity,
                                                 a = lower,
                                                 b = upper)
    
    state$demography$transition_matrix <- transition_matrix
    
    state
    
  }
  
  demographic_dynamics
  
}

environmental_stochasticity_list <- function (transition_matrix_list) {

  demographic_dynamics <- function (state, timestep) {
    
    state$demography$transition_matrix <- transition_matrix_list[[timestep]]
    
    state
    
  }
  
  demographic_dynamics
  
}


# the main exported user function - to run an experiment
experiment <- function (habitat,
                        demography,
                        dynamics = simple_dynamics,
                        timesteps = 100) {
  
  
  # build the state object
  state <- list(habitat = habitat,
                demography = demography)
  
  # iterate the system
  state <- iterate_system(state, seq_len(timesteps), dynamics)
  
  state
  
}

# user friendly functions to construct the demography and habitate objects
demography <- function (transition_matrix) {
  object <- list(transition_matrix = transition_matrix)
  object <- set_class(object, "demography")
  object
}

habitat <- function (locations, initial_population, carrying_capacity = NULL) {
  object <- list(population = initial_population,
                 carrying_capacity = carrying_capacity,
                 locations = locations)
  object <- set_class(object, "habitat")
  object
}

build_dynamics <- function (transition,
                            dispersal,
                            demographic_dynamics,
                            habitat_dynamics,
                            order = c("habitat_dynamics",
                                      "demographic_dynamics",
                                      "transition",
                                      "dispersal")) {
  # get all the functions in a list, in the required order
  check_transition_order(order)
  functions <- lapply(order, get, envir = environment())
  functions
}


simple_dynamics <- build_dynamics(transition = simple_transition,
                                  dispersal = simple_dispersal(1),
                                  demographic_dynamics = no_demographic_dynamics,
                                  habitat_dynamics = no_habitat_dynamics,
                                  order = c("habitat_dynamics",
                                            "demographic_dynamics",
                                            "transition",
                                            "dispersal")) 


# ~~~~~~~~~~~~
# test run

set.seed(1)
n_patches <- 100
n_stages <- 3

# carrying capacity of *total* population
capacity <- rpois(n_patches, 5)
pop <- fake_population(capacity, n_stages)
locs <- matrix(runif(n_patches * 2), ncol = 2)
trans <- fake_transition_matrix(n_stages)


# create the two objects defining the system (no carrying capacities)
hab <- habitat(locs, pop)
demog <- demography(trans)

# simple experiment (no effect of carrying capacity)
out <- experiment(hab, demog)


# ~~~~
# add carrying capacity

my_dynamic <- build_dynamics(transition = simple_transition_capped,
                          dispersal = simple_dispersal(1),
                          habitat_dynamics = no_habitat_dynamics,
                          demographic_dynamics = no_demographic_dynamics)

out <- experiment(hab, demog, my_dynamic)
# ut oh, we forgot to supply a carrying capacity!

hab <- habitat(locs, pop, capacity)
out <- experiment(hab, demog, my_dynamic)
pops <- sapply(out,
               function(state) state$habitat$population[1, ])

plot(pops[1, ], type = "l", ylim = c(0, 3))
lines(pops[2, ], col = "red")
lines(pops[3, ], col = "blue")
# that's better

# ~~~~
# every 5 years, the population is culled by a half
my_transition <- function (state, timestep) {
  state <- simple_transition_capped(state, timestep)
  
  if (timestep %% 5 == 0) {
    state$habitat$population <- state$habitat$population * 0.5
  }
  
  state
  
}


my_dynamic <- build_dynamics(transition = my_transition,
                             dispersal = simple_dispersal(1),
                             habitat_dynamics = no_habitat_dynamics,
                             demographic_dynamics = no_demographic_dynamics)

out <- experiment(hab, demog, my_dynamic)

pops <- vapply(out,
               function(state) state$habitat$population[1, 3],
               1)

plot(pops, type = "l")


# ~~~~
# add environmental stochasticity, in a boring way
out <- experiment(hab, demog,
                  transition = simple_transition_capped,
                  demographic_dynamics = environmental_stochasticity(trans, 0.5))

# ~~~~
# add environmental stochasticity, from a predefined list of tranistion matrices
transitions <- replicate(100,
                         fake_transition_matrix(n_stages),
                         simplify = FALSE)

out <- experiment(hab, demog,
                  transition = simple_transition_capped,
                  demographic_dynamics = environmental_stochasticity_list(transitions))

# out$demography$transition_matrix
# transitions[[100]]
