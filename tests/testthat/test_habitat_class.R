context('habitat-class')

test_that('habitat classes work', {
  library(raster)

  expect_true(inherits(habitat(egk_hab, egk_k),"habitat"))
  
  expect_true(is.habitat(habitat(egk_hab, egk_k)))

  print(habitat(egk_hab, egk_k))

})
