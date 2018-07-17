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
  r <- raster(vals=1, nrows=150, ncols=150, res=c(5,5), crs=('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
  
  hab.suit <- r*sample(seq(0,1,.01), ncell(r), replace=TRUE)
  
  r2 <- r
  r2[] <- 0
  cells <- sample(c(1:ncell(r2)), 10)
  r2[c(adjacent(hab.suit, cells, directions=8, pairs=FALSE),cells)]  <- 10
  r3 <- r2#*hab.suit
  
  pop <- stack(r3*2,r3*2,r3*2,r3*2)
  
  hab.k <- hab.suit*10
  
  disp.bar <- hab.suit*0
  disp.bar[cellFromCol(disp.bar,ncol(disp.bar)/2)] <- 1
  disp.bar2 <- hab.suit*0
  disp.bar2[sampleRandom(disp.bar2, size=700, na.rm=TRUE, sp=TRUE)] <- 1
  
  params <- list(
    dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=1,'Stage_2-3'=1,'Stage_3+'=0),
    dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
    dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0)#,
    #barriers_map=disp.bar,
    #use_barriers=TRUE
  )
  
  params2 <- list(
    dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
    dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
    dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
    barriers_map=disp.bar,
    use_barriers=TRUE,
    barrier_type=1
  )
  
  params3 <- list(
    dispersal_distance=list('Stage_0-1'=0,'Stage_1-2'=10,'Stage_2-3'=10,'Stage_3+'=0),
    dispersal_kernel=list('Stage_0-1'=0,'Stage_1-2'=exp(-c(0:9)^1/3.36),'Stage_2-3'=exp(-c(0:9)^1/3.36),'Stage_3+'=0),
    dispersal_proportion=list('Stage_0-1'=0,'Stage_1-2'=0.35,'Stage_2-3'=0.35*0.714,'Stage_3+'=0),
    barriers_map=disp.bar2,
    use_barriers=TRUE,
    barrier_type=1
  )
  
  expect_true(inherits(build_demography(mat,dispersal_parameters=params3),"demography"))
  
  plot(build_demography(mat,dispersal_parameters=params3))
  
  expect_error(build_demography(mat[c(1:2),c(1:3)],dispersal_parameters=params3))
  
  expect_true(inherits(build_demography(mat,dispersal_parameters=params3),"demography"))
  
  expect_true(inherits(build_demography(transition_matrix=mat,
                                scale='local',
                                habitat_suitability=hab.suit,
                                dispersal_parameters=params),
                       "demography")
              )
  
  expect_error(build_demography(transition_matrix=mat,
                                        scale='local',
                                        dispersal_parameters=params)
               )
  
  mat2 <- mat
  mat2[1,2] <- NA
  expect_error(build_demography(mat2,dispersal_parameters=params3))
  
  expect_error(build_demography(as.vector(mat),dispersal_parameters=params3))
  
  print(build_demography(mat,dispersal_parameters=params3))
  
  summary(build_demography(mat,dispersal_parameters=params3))
  
  plot(build_demography(mat,dispersal_parameters=params3))
  
  expect_true(is.demography(build_demography(mat,dispersal_parameters=params3)))
  
  #expect_error(as.demography(c(1,2,3)))

})
 