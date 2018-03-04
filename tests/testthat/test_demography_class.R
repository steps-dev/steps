context('demography-class')

test_that('demography classes work', {
  library(raster)
  library(rgdal)
    
  # the types of demography
  mat <- matrix(c(0.000,0.000,0.302,0.302,
                  0.940,0.000,0.000,0.000,
                  0.000,0.884,0.000,0.000,
                  0.000,0.000,0.793,0.793),
                nrow = 4, ncol = 4, byrow = TRUE)
  colnames(mat) <- rownames(mat) <- c('Stage_0-1','Stage_1-2','Stage_2-3','Stage_3+')
    
  # the types of demography
  mat_sd <- matrix(c(0.000,0.000,1,1,
                     1,0.000,0.000,0.000,
                     0.000,1,0.000,0.000,
                     0.000,0.000,1,1),
                   nrow = 4, ncol = 4, byrow = TRUE)
  colnames(mat_sd) <- rownames(mat_sd) <- c('Stage_0-1','Stage_1-2','Stage_2-3','Stage_3+')

  # the types of habitat attributes
  r <- raster(vals=1, nrows=10, ncols=10, res=100, crs=('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  r2 <- r
  r2[] <- round(runif(ncell(r2),0,20),0)
  
  r3 <- projectRaster(r2, crs=("+init=epsg:4326")) 
  
  ss.dist <- c(0.3,0.4,0.3,0.5)
  
  hab.suit <- r
  
  hab.suit.s <- stack(r,r,r,r,r)

  hab.pop <- r2
  
  hab.pop.s <- stack(r2,r2,r2,r2,r2)
  
  hab.pop.n <- 4
  hab.pop.n2 <- c(4,4,4,4)
  
  hab.pop.sp <- sampleRandom(r, size=3, na.rm=TRUE, sp=TRUE) 
  hab.pop.sp@data <- as.data.frame(t(rmultinom(3, size = 100, prob = c(0.8,0.2,0.1))))
  hab.pop.sp2 <- spTransform(hab.pop.sp, CRS("+init=epsg:28355"))
  
  hab.k <- r2*2
  
  hab.k.func <- function(x) x*0.8
  
  habitat <- as.habitat(list(as.habitat_suitability(hab.suit), as.populations(hab.pop.n2), as.carrying_capacity(hab.k)))
  habitat_popNA <- as.habitat(list(as.habitat_suitability(hab.suit), as.populations(hab.pop.n2), as.carrying_capacity(hab.k)))
  habitat_altK <- as.habitat(list(as.habitat_suitability(hab.suit), as.populations(hab.pop.n2), as.carrying_capacity(hab.k*.01)))
  populations(habitat_popNA)[[2]][c(1,5)] <- NA
    
  # check as.demography won't handle a silly function
  expect_error(as.demography(mat[,-1]))
  expect_error(as.demography(mat[-1,]))
  expect_error(as.demography(mat[-1,], type='local', hab.suit))
  
  # check they have the right class
  expect_s3_class(as.demography(mat), 'demography')
  expect_s3_class(as.demography(mat, type='local', habsuit=hab.suit), 'demography')
  expect_s3_class(as.demography(mat, type='local', habsuit=hab.suit.s), 'demography')
  
  # check for proper require inputs
  expect_error(as.demography(mat, type='local'))
  expect_error(as.demography(mat, type='local', 1))

  # check is.demography works on demographies
  expect_true(is.demography(as.demography(mat)))
    
  # check for NA values whilst estimating demography
  expect_error(estimate_demography(as.demography(mat),habitat_popNA, time_step=1))
  
  # check output of estimate demography
  expect_true(inherits(estimate_demography(as.demography(mat),habitat, time_step=1), 'list'))
  expect_true(inherits(estimate_demography(as.demography(mat),habitat, time_step=1)[[1]], 'RasterLayer'))
  
  expect_true(inherits(estimate_demography(as.demography(mat),habitat_altK, time_step=1), 'list'))
  expect_true(inherits(estimate_demography(as.demography(mat),habitat_altK, time_step=1)[[1]], 'RasterLayer'))
  
  expect_true(inherits(estimate_demography(as.demography(mat),habitat, time_step=1, seed=123), 'list'))
  expect_true(inherits(estimate_demography(as.demography(mat),habitat, time_step=1, seed=123)[[1]], 'RasterLayer'))

  expect_true(inherits(estimate_demography(as.demography(mat, type="local", hab.suit), habitat, time_step=1), 'list'))
  expect_true(inherits(estimate_demography(as.demography(mat, type="local", hab.suit), habitat, time_step=1)[[1]], 'RasterLayer'))
  
  expect_true(inherits(estimate_demography(as.demography(mat, type="local", hab.suit),habitat, time_step=1, seed=123), 'list'))
  expect_true(inherits(estimate_demography(as.demography(mat, type="local", hab.suit),habitat, time_step=1, seed=123)[[1]], 'RasterLayer'))
      
  
  summary.demography(as.demography(mat))
  plot.demography(as.demography(mat))
})
  