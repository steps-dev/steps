context('demography-class')

test_that('demography classes work', {
  library(raster)
  library(rgdal)

  expect_true(inherits(demography(egk_mat),"demography"))
  
  plot(demography(egk_mat))
  
  expect_error(demography(egk_mat[c(1:2),c(1:3)]))
  
  expect_true(inherits(demography(transition_matrix=egk_mat,
                                scale='local',
                                habitat_suitability=egk_hab),
                       "demography")
              )
  
  expect_error(demography(transition_matrix=egk_mat,
                                        scale='local')
               )
  
  mat2 <- egk_mat
  mat2[1,2] <- NA
  expect_error(demography(mat2))
  
  expect_error(demography(as.vector(egk_mat)))
  
  print(demography(egk_mat))
  
  summary(demography(egK_mat))
  
  plot(demography(egk_mat))
  
  expect_true(is.demography(demography(egk_mat)))

})
 