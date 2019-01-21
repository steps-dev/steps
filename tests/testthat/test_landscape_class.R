context('landscape-class')

test_that('landscape classes work', {
  library(raster)
  
  expect_true(inherits(landscape(population = egk_pop,
                                 suitability = egk_hab,
                                 carrying_capacity = egk_k),
                       "landscape"))
  
  expect_true(is.landscape(landscape(population = egk_pop,
                                     suitability = egk_hab,
                                     carrying_capacity = egk_k)))
  
  expect_error(landscape(population = NULL,
                         suitability = egk_hab,
                         carrying_capacity = egk_k))
  
  pop2 <- egk_pop
  res(pop2) <- 10
  expect_error(landscape(population = pop2,
                         suitability = egk_hab,
                         carrying_capacity = egk_k))
  
  pop3 <- egk_pop
  extent(pop3) <- c(-160, 180, -90, 90)
  res(pop3) <- 5
  expect_error(landscape(population = pop3,
                         suitability = egk_hab,
                         carrying_capacity = egk_k))
  
  print(landscape(population = egk_pop,
                  suitability = egk_hab,
                  carrying_capacity = egk_k))
  
})
 