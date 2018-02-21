context('habitat-class')

test_that('habitat classes work', {
  
  # the types of habitat attributes
  r <- raster(vals=1, nrows=10, ncols=10, res=100, crs=('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  hab.suit <- r
  hab.suit.s <- stack(r,r,r,r,r)
  
  r2 <- r
  r2[] <- round(runif(ncell(r2),0,20),0)
  hab.pop <- r2
  hab.pop.s <- stack(r2,r2,r2,r2,r2)
  hab.pop.n <- 4
  
  hab.k <- r2*2
  hab.k.func <- function(x) x*0.8
  
  expect_true(inherits(as.habitat_suitability(hab.suit),"habitat_suitability"))
  expect_true(inherits(as.habitat_suitability(hab.suit.s),"habitat_suitability"))

  expect_error(as.habitat_suitability(1))

  expect_true(inherits(as.populations(hab.pop),"populations"))
  expect_true(inherits(as.populations(hab.pop.s),"populations"))
  expect_true(inherits(as.populations(hab.pop.n),"populations"))
  
  expect_error(as.populations("a"))
  
  expect_true(inherits(as.carrying_capacity(hab.k),"carrying_capacity"))
  expect_true(inherits(as.carrying_capacity(hab.k.func),"carrying_capacity"))
  
  expect_error(as.carrying_capacity(1))
  
  expect_error(as.habitat(c(1,2,3)))
  expect_error(as.habitat(list(1,2)))
  expect_error(as.habitat(list(as.populations(hab.pop),as.habitat_suitability(hab.suit),as.carrying_capacity(hab.k))))
  expect_error(as.habitat(list(as.habitat_suitability(hab.suit),as.carrying_capacity(hab.k),as.populations(hab.pop))))
  expect_error(as.habitat(list(as.habitat_suitability(hab.suit),as.carrying_capacity(hab.k),as.carrying_capacity(hab.k))))

  expect_error(list2habitat())
    
  # check they have the right class
  expect_s3_class(as.habitat(list(as.habitat_suitability(hab.suit),as.populations(hab.pop.n),as.carrying_capacity(hab.k)), 'habitat')
  
  # check is.habitat works on habitats
  expect_true(is.habitat(as.habitat(list(as.habitat_suitability(hab.suit),as.populations(hab.pop.n),as.carrying_capacity(hab.k))))
  
})