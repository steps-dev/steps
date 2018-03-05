context('dispersal-class')

test_that('dispersal functions work', {
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
  r <- raster(vals=1, nrows=10, ncols=10, crs=('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  r2 <- r
  r2[] <- round(runif(ncell(r2),0,20),0)
  
  r3 <- projectRaster(r2, crs=("+init=epsg:4326")) 
  
  ss.dist <- c(0.3,0.4,0.3,0.5)
  
  hab.suit <- r
  
  hab.suit.s <- stack(r,r,r,r,r)

  hab.pop <- r2
  
  hab.pop.s <- stack(r2,r2,r2,r2)
  
  hab.pop.n <- 4
  hab.pop.n2 <- c(4,4,4,4)
  
  hab.pop.sp <- sampleRandom(r, size=3, na.rm=TRUE, sp=TRUE) 
  hab.pop.sp@data <- as.data.frame(t(rmultinom(3, size = 100, prob = c(0.8,0.2,0.1))))
  hab.pop.sp2 <- spTransform(hab.pop.sp, CRS("+init=epsg:28355"))
  
  hab.k <- r2*2
  
  hab.k.func <- function(x) x*0.8
  
  habitat <- as.habitat(list(as.habitat_suitability(hab.suit), as.populations(hab.pop.n2), as.carrying_capacity(hab.k)))
  habitat2 <- as.habitat(list(as.habitat_suitability(hab.suit.s), as.populations(hab.pop.n2), as.carrying_capacity(hab.k)))
  habitat3 <- as.habitat(list(as.habitat_suitability(hab.suit.s), as.populations(hab.pop.n2), as.carrying_capacity(hab.k*.01)))

  disp.bar <- hab.suit*0
  disp.bar[cellFromCol(disp.bar,ncol(disp.bar)/2)] <- 1
  disp.bar2 <- hab.suit*0
  disp.bar2[sampleRandom(disp.bar2, size=25, na.rm=TRUE, sp=TRUE)] <- 1
  
  disp.bar.s <- stack(mget(rep("disp.bar",3)))
  
  params <- as.dispersal(
    list(
      dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
      dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
      dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0)
    )
  )
  
  paramsBT <- as.dispersal(
    list(
      dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
      dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
      dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
      barrier_type=1
    )
  )
  
  paramsDS <- as.dispersal(
    list(
      dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
      dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
      dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
      dispersal_steps=2
    )
  )

  paramsUB <- as.dispersal(
    list(
      dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
      dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
      dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
      use_barriers=TRUE
    )
  )

  paramsBM <- as.dispersal(
    list(
      dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
      dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
      dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
      barriers_map=disp.bar
    )
  )

  paramsBM2 <- as.dispersal(
    list(
      dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
      dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
      dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
      barriers_map=disp.bar2
    )
  )
    
  paramsBMS <- as.dispersal(
    list(
      dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
      dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
      dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
      barriers_map=disp.bar.s
    )
  )
          
  # check for proper require inputs
  expect_error(as.dispersal(1))
  expect_error(as.dispersal(list(1)))
  expect_error(as.dispersal(list('a'=1,'b'=2,'c'=3)))
  
  # check they have the right class
  expect_s3_class(
    as.dispersal(
      list(
        dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
        dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
        dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0)
      )
    ),
    'dispersal')

  # check is.dispersal works on dispersal objects
  expect_true(
    is.dispersal(
      as.dispersal(
        list(
          dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
          dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
          dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0)
        )
      )
    )
  )
    
  # check missing method type in dispersal call
  expect_error(
    dispersal(params=params, habitat=habitat) 
  )
  
  expect_error(
    dispersal(params=list('a'=1,'b'=2,'c'=3), habitat=habitat, method='ca') 
  )

  expect_error(
    dispersal(params=params, habitat=hab.suit, method='ca') 
  )
  

    
  # check output of dispersal - cellular automata
  expect_true(inherits(dispersal(params, habitat, method='ca', time_step=1), 'list'))

  expect_true(inherits(dispersal_core_ca(params, habitat, time_step=1), 'list'))
  expect_true(inherits(dispersal_core_ca(params, habitat, time_step=1)[[1]], 'RasterLayer'))
  
  expect_true(inherits(dispersal_core_ca(paramsBT, habitat, time_step=1), 'list'))
  expect_true(inherits(dispersal_core_ca(paramsDS, habitat, time_step=1), 'list')) 
  expect_true(inherits(dispersal_core_ca(paramsUB, habitat, time_step=1), 'list'))
  
  expect_true(inherits(dispersal_core_ca(params, habitat2, time_step=1), 'list'))
  
  expect_true(inherits(dispersal_core_ca(params, habitat3, time_step=1), 'list'))

  expect_true(inherits(dispersal_core_ca(paramsBM, habitat, time_step=1), 'list'))
  expect_true(inherits(dispersal_core_ca(paramsBM2, habitat, time_step=1), 'list'))
  expect_true(inherits(dispersal_core_ca(paramsBMS, habitat, time_step=1), 'list'))

  # check output of dispersal - fast fourier transformation 
  expect_true(inherits(dispersal(params, habitat, method='fft'), 'list'))
  expect_true(inherits(dispersal_core_fft(params, habitat), 'list'))
  expect_true(inherits(dispersal_core_fft(params, habitat)[[1]], 'RasterLayer'))  
  
  print(
    as.dispersal(
      list(
        dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
        dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
        dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0)
      )
    )
  )

})
  