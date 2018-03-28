## ---- message = FALSE----------------------------------------------------
library(steps)
library(raster)
library(future)
library(rasterVis)

koala.trans.mat <- matrix(c(0.000,0.000,0.302,0.302,
                              0.940,0.000,0.000,0.000,
                              0.000,0.884,0.000,0.000,
                              0.000,0.000,0.793,0.793),
                            nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat) <- rownames(koala.trans.mat) <- c('Juveniles','Sub_Adults','Adults','Super_Adults')

koala.trans.mat.es <- matrix(c(0.000,0.000,1,1,
                                 1,0.000,0.000,0.000,
                                 0.000,1,0.000,0.000,
                                 0.000,0.000,1,1),
                               nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat.es) <- rownames(koala.trans.mat.es) <- c('Juveniles','Sub_Adults','Adults','Super_Adults')


## ---- message = FALSE----------------------------------------------------
library(raster)

# read in spatial habitat suitability raster
koala.hab.suit <- raster(system.file("extdata","Koala_HabSuit.tif", package="steps"))

# rescale habitat raster to relative suitability (or likelihood of occurrence) between 0 and 1.
koala.hab.suit <- (koala.hab.suit - cellStats(koala.hab.suit, min)) / (cellStats(koala.hab.suit, max) - cellStats(koala.hab.suit, min))

# rename raster so that it is tracked throughout the experiment and simulation
names(koala.hab.suit) <- "Habitat"

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(koala.hab.suit, box = FALSE, axes = FALSE)

## ---- message = FALSE----------------------------------------------------

# create carrying capacity layer using the habitat suitability raster (or provide a custom one)
koala.hab.k <- ceiling(koala.hab.suit * 60)

# rename raster so that it is tracked throughout the experiment and simulation
names(koala.hab.k) <- "Carrying Capacity"

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(koala.hab.k, box = FALSE, axes = FALSE)

## ---- message = FALSE----------------------------------------------------

# create population layers using the carrying capacity raster (or provide a custom ones)
koala.pop <- stack(replicate(4, ceiling(koala.hab.k * 0.2)))

# rename stack so that the layers are tracked throughout the experiment and simulation
names(koala.pop) <- colnames(koala.trans.mat)

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spplot(koala.pop)

## ---- message = FALSE----------------------------------------------------

# define all of the dispersal distances
dispersal_distance <- list('Juveniles'=0,'Sub_Adults'=10,'Adults'=10,'Super_Adults'=0)

# define the dispersal kernals for each life-stage
dispersal_kernel <- list('Juveniles'=0,'Sub_Adults'=exp(-c(0:9)^1/3.36),'Adults'=exp(-c(0:9)^1/3.36),'Super_Adults'=0)

# define all of the dispersal proportions
dispersal_proportion <- list('Juveniles'=0,'Sub_Adults'=0.35,'Adults'=0.35*0.714,'Super_Adults'=0)

# combine all of the parameters in a list
koala.disp.param <- list(dispersal_distance=dispersal_distance,
                                      dispersal_kernel=dispersal_kernel,
                                      dispersal_proportion=dispersal_proportion
                         )
