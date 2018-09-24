context('population-class')

test_that('population classes work', {
  library(raster)

  expect_true(inherits(population(egk_pop),"population"))
  
  expect_true(is.population(population(egk_pop)))

  print(population(egk_pop))
  
})
