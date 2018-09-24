context('state-class')

test_that('state classes work', {
  library(raster)
  
  b_pop <- population(egk_pop)
  
  b_hab <- habitat(habitat_suitability = egk_hab,
                       carrying_capacity = egk_k)

  b_dem <- demography(transition_matrix = egk_mat)
  
  expect_true(inherits(state(population = b_pop,
                             habitat = b_hab,
                             demography = b_dem),
                       "state")
              )
  
  expect_true(is.state(state(population = b_pop,
                             habitat = b_hab,
                             demography = b_dem)
                          )
              )
  
  expect_error(build_state(population = b_pop,
                           habitat = b_hab)
               )
 
  pop2 <- egk_pop
  res(pop2) <- 10
  b_pop2 <- population(pop2)
  expect_error(state(population = b_pop2,
                     habitat = b_hab,
                     demography = b_dem)
               )
  
  pop3 <- egk_pop
  extent(pop3) <- c(-160, 180, -90, 90)
  res(pop3) <- 5
  b_pop3 <- population(pop3)
  expect_error(state(population = b_pop3,
                           habitat = b_hab,
                           demography = b_dem)
               )

  mat2 <- egk_mat[c(1:2),c(1:2)]
  b_dem2 <- build_demography(transition_matrix = mat2)
  expect_error(state(population = b_pop,
                     habitat = b_hab,
                     demography = b_dem2)
               )
    
  
  print(state(population = b_pop,
              habitat = b_hab,
              demography = b_dem)
        )

})
 