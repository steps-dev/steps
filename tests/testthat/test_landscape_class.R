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

  egk_pop2 <- egk_pop
  res(egk_pop2) <- 250
  expect_error(landscape(population = egk_pop2,
                         suitability = egk_hab,
                         carrying_capacity = NULL))
  
  egk_pop3 <- egk_pop
  extent(egk_pop3) <- extent(329000, 337000, 5817000, 5834500)
  res(egk_pop3) <- 500
  expect_error(landscape(population = egk_pop3,
                         suitability = egk_hab,
                         carrying_capacity = NULL))
  
  egk_pop4 <- egk_pop
  egk_pop4[c(1:20)] <- NA
  expect_error(landscape(population = egk_pop4,
                         suitability = egk_hab,
                         carrying_capacity = NULL))
  
  # print(landscape(population = egk_pop,
  #                 suitability = egk_hab,
  #                 carrying_capacity = egk_k))
  
})
 