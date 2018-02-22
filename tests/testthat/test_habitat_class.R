context('habitat-class')

test_that('habitat classes work', {
  library(raster)
  library(rgdal)

  # the types of habitat attributes
  r <- raster(vals=1, nrows=10, ncols=10, res=100, crs=('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  hab.suit <- r
  hab.suit.s <- stack(r,r,r,r,r)
  
  r2 <- r
  r2[] <- round(runif(ncell(r2),0,20),0)
  hab.pop <- r2
  hab.pop.s <- stack(r2,r2,r2,r2,r2)
  hab.pop.n <- 4
  
  hab.pop.sp <- sampleRandom(r, size=3, na.rm=TRUE, sp=TRUE) 
  hab.pop.sp@data <- as.data.frame(t(rmultinom(3, size = 100, prob = c(0.8,0.2,0.1))))
  hab.pop.sp2 <- spTransform(hab.pop.sp, CRS("+init=epsg:28355"))
  
  hab.k <- r2*2
  hab.k.func <- function(x) x*0.8
  
  r3 <- projectRaster(r2, crs=("+init=epsg:4326"))
  
  hab.suit.s.NA <- hab.suit.s
  hab.suit.s.NA[[1]][1,c(1,2)] <- "a"
  
  hab.suit.s.0 <- hab.suit.s
  hab.suit.s.0[[1]][1,c(1,2)] <- 0
  
  expect_identical(attr(as.habitat_suitability(hab.suit), "habitat"), "habitat_suitability")
  expect_identical(attr(as.habitat_suitability(hab.suit.s), "habitat"), "habitat_suitability")

  expect_error(as.habitat_suitability(1))

  expect_error(habitat_suitability(1))
  expect_true(inherits(habitat_suitability(as.habitat(list(as.habitat_suitability(hab.suit),as.populations(hab.pop.n),as.carrying_capacity(hab.k)))),c("RasterLayer")))
  
  expect_identical(attr(as.populations(hab.pop), "habitat"), "populations")
  expect_identical(attr(as.populations(hab.pop.s), "habitat"), "populations")
  expect_identical(attr(as.populations(hab.pop.n), "habitat"), "populations")
  expect_identical(attr(as.populations(hab.pop.sp), "habitat"), "populations")

  expect_error(as.populations("a"))

  expect_error(populations(1))
  expect_true(inherits(populations(as.habitat(list(as.habitat_suitability(hab.suit),as.populations(hab.pop.n),as.carrying_capacity(hab.k)))),c("list")))

  #expect_error(populations(as.habitat(list(as.habitat_suitability(hab.suit),as.populations(hab.pop.n),as.carrying_capacity(hab.k)))) <- 1)
   
  expect_error(populations2rasterbrick(hab.pop.sp2,hab.suit))
  expect_true(inherits(populations2rasterbrick(hab.pop.sp, hab.suit),c("RasterBrick")))

  expect_identical(attr(as.carrying_capacity(hab.k), "habitat"), "carrying_capacity")
  expect_identical(attr(as.carrying_capacity(hab.k.func), "habitat"), "carrying_capacity")

  expect_error(as.carrying_capacity(1))

  expect_error(carrying_capacity(1))
  expect_true(inherits(carrying_capacity(as.habitat(list(as.habitat_suitability(hab.suit),as.populations(hab.pop.n),as.carrying_capacity(hab.k)))),c("RasterLayer")))
  
  expect_error(list2habitat())
  
  expect_error(area_check(hab.suit.s.NA))
  expect_error(area_check(hab.suit.s.0))
  
  expect_true(inherits(area_of_region(r3),c("RasterLayer")))

  expect_error(as.habitat(c(1,2,3)))
  expect_error(as.habitat(list(1,2)))
  expect_error(as.habitat(list(as.populations(hab.pop),as.habitat_suitability(hab.suit),as.carrying_capacity(hab.k))))
  expect_error(as.habitat(list(as.habitat_suitability(hab.suit),as.carrying_capacity(hab.k),as.populations(hab.pop))))
  expect_error(as.habitat(list(as.habitat_suitability(hab.suit),as.carrying_capacity(hab.k),as.carrying_capacity(hab.k))))

  # check they have the right class
  expect_s3_class(as.habitat(list(as.habitat_suitability(hab.suit),as.populations(hab.pop.n),as.carrying_capacity(hab.k))), 'habitat')
  expect_s3_class(as.habitat(list(as.habitat_suitability(hab.suit.s),as.populations(hab.pop.n),as.carrying_capacity(hab.k))), 'habitat')
  
  # check is.habitat works on habitats
  expect_true(is.habitat(as.habitat(list(as.habitat_suitability(hab.suit),as.populations(hab.pop.n),as.carrying_capacity(hab.k)))))

  print(as.habitat(list(as.habitat_suitability(hab.suit),as.populations(hab.pop.n),as.carrying_capacity(hab.k))))
  
})
  