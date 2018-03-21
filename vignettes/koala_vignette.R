## ------------------------------------------------------------------------
koala.trans.mat <- matrix(c(0.000,0.000,0.302,0.302,
                              0.940,0.000,0.000,0.000,
                              0.000,0.884,0.000,0.000,
                              0.000,0.000,0.793,0.793),
                            nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat) <- rownames(koala.trans.mat) <- c('Juveniles','Sub-Adults','Adults','Super-Adults')

koala.trans.mat.es <- matrix(c(0.000,0.000,1,1,
                                 1,0.000,0.000,0.000,
                                 0.000,1,0.000,0.000,
                                 0.000,0.000,1,1),
                               nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat.es) <- rownames(koala.trans.mat.es) <- c('Juveniles','Sub-Adults','Adults','Super-Adults')


## ------------------------------------------------------------------------
library(raster)
# read in spatial habitat suitability raster
koala.hab.suit <- raster(system.file("extdata","Koala_HabSuit.tif", package="steps"))

# rescale habitat raster to relative suitability (or likelihood of occurrence) between 0 and 1.
koala.hab.suit <- (koala.hab.suit - cellStats(koala.hab.suit, min)) / (cellStats(koala.hab.suit, max) - cellStats(koala.hab.suit, min))

# rename raster so that it is traked throughout the experiment and simulation
names(koala.hab.suit) <- "Habitat"
plot(koala.hab.suit, box = FALSE, axes = FALSE)

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

