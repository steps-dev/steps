## test the dispersal c++ function.
library(raster)
set.seed(42) #grandstand
xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
Dd <- as.matrix(dist(xy))
w <- exp(-1/nrow(xy) * Dd)
Ww <- chol(w)
xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
coordinates(xy) <- ~x+y
r <- rasterize(xy, raster(points2grid(xy)), 'z')
proj4string(r) <- '+init=epsg:4283'
r <- disaggregate(r,10,method='bilinear')

xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
Dd <- as.matrix(dist(xy))
w <- exp(-1/nrow(xy) * Dd)
Ww <- chol(w)
xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
coordinates(xy) <- ~x+y
r2 <- rasterize(xy, raster(points2grid(xy)), 'z')
proj4string(r2) <- '+init=epsg:4283'
r2 <- disaggregate(r2,10,method='bilinear')

params <- list(dispersal_kernal = list('larvae'=exp(-c(0:4)), 'juvenile'=0, 'adult'=exp(-c(0:4)*3)),
               dispersal_proportion = list('larvae'=0.2, 'juvenile'=0,'adult'=0.6))  

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
r[] <- range01(r[])
r2[] <- range01(r2[])

starting_population_state <- r*r2*100
starting_population_state[starting_population_state[]<45]<-0

potiential_carrying_capacity <- r*150
# potiential_carrying_capacity<-raster::flip(potiential_carrying_capacity,direction = 'y')
habitat_suitability_map <- r
# habitat_suitability_map<-raster::flip(habitat_suitability_map,direction = 'y')
barriers_map <- r
# barriers_map[runif(51*51)>.95] <- NA
barriers_map[!is.na(barriers_map)] <-0
barriers_map[250:255,125:375]<-1
# barriers_map[,26:51]<-0

barrier_type = 0
use_barrier = TRUE
dispersal_steps = 1
dispersal_distance = 10
dispersal_kernel = exp(-c(0:9)*.1)
dispersal_proportion = .3

dispersal_list <- list()

# starting_population_state[starting_population_state[]<40]<-0
dispersal_list[[1]] <- as.matrix(starting_population_state)
for (i in 1:11){
dispersal_list[[i+1]] <- dhmpr::rcpp_dispersal(as.matrix(dispersal_list[[i]]), as.matrix(potiential_carrying_capacity), as.matrix(habitat_suitability_map), as.matrix(barriers_map), barrier_type, use_barrier, dispersal_steps, dispersal_distance, dispersal_kernel, dispersal_proportion)[[1]]
print(i)
}

rasters <- lapply(dispersal_list,raster)
rb <- brick(rasters)
library(viridis)
animation::saveGIF(animate(rb,pause=.001,main='dispersal',n=1,col=rev(viridis::viridis(100))), movie.name = "test_dispersal_t50_500x500_w_barrier.gif",loop=TRUE)
