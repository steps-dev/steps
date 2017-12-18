## message for Casey: this is some code to see if we can set up the habitats, populations and other functions as sparseMatricies.
## hopefully if we can, we can go about converting internal R functions to run on sparseMatrices and away we go.

library(Matrix)
library(raster)
library(dhmpr)
set.seed(42)
xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
Dd <- as.matrix(dist(xy))
w <- exp(-1/nrow(xy) * Dd)
Ww <- chol(w)
xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
coordinates(xy) <- ~x+y
r <- rasterize(xy, raster(points2grid(xy)), 'z')
proj4string(r) <- '+init=epsg:4283'
r[] <- scales::rescale(r[],to=c(0,1))

## create a habitat from a list containing a habitat suitability raster and numeric values for population and carrying capacity.
hsm <- as.habitat_suitability(r)
pops <- as.populations(c(80,20,10))
cc <- as.carrying_capacity(300)

features <- list(hsm,pops,cc)
habitat <- as.habitat(features)

## check out size of different objects
sM<-Matrix(getValues(habitat_suitability(habitat)))
object.size(sM)
object.size(habitat_suitability(habitat))


sM<-Matrix(getValues(stack(populations(habitat))))
object.size(sM)
object.size(populations(habitat))

