#' @title habitat objects
#' @name habitat
#' @rdname habitat
#' @description Underlying habitat for dlmpr.

#' @rdname habitat
#' @name patchify
#' @param x a binary Raster layer (0 or NA for background, and 1 for areas to be clumped)
#' @param distance the neighbourhood distance. Patches that occur within this distance of
#'  one another will be clumped. This should be in the units of the CRS given
#'  in p4s. If this is zero, the function will identify patches defined by 
#'  8-connectivity (i.e. queen's case).
#' @param p4s an equal-area projection appropriate for the region. This will be used when 
#'      calculating inter-patch distances. The returned objects will be in the 
#'      original coordinate system.
#' @param givedist should the distance matrix be returned? (logical). Distances are in the 
#'           units of p4s, and are shortest cartesian distances between patches.
#' @author John Baumgartner
#' @export

#' @importFrom raster raster
#' @importFrom sp CRS
#' @importFrom methods is
#' @importFrom sp spDists
#' @importFrom raster xyFromCell
#' @importFrom raster clump
#' @importFrom raster getValues
#' @importFrom raster ncell
#' @importFrom raster freq
#' @importClassesFrom raster Raster
#' @importClassesFrom sp CRS

#' @examples
#' library(raster)
#' set.seed(42)
#' xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
#' Dd <- as.matrix(dist(xy))
#' w <- exp(-1/nrow(xy) * Dd)
#' Ww <- chol(w)
#' xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
#' coordinates(xy) <- ~x+y
#' r <- rasterize(xy, raster(points2grid(xy)), 'z')
#' r2 <- raster(r)
#' res(r2) <- 0.01
#' r2 <- resample(r, r2)
#' proj4string(r2) <- '+init=epsg:4283'
#' rthr <- r2 > quantile(r2, 0.9)
#' foo <- patchify(rthr, distance=1000, '+init=epsg:3577')
#' plot(foo$patchrast)
#' text(foo$coords[, 2:3], labels = foo$coords[, 1])

patchify <- function(x, distance, p4s, givedist=TRUE) {
  if(!is(x, 'Raster')) x <- raster::raster(x)
  if(!is(p4s, 'CRS')) p4s <- sp::CRS(p4s)
  if(base::is.na(sp::proj4string(x))) stop(base::substitute(x), ' lacks a CRS.')
  rc <- raster::clump(x,directions = 4)
  clump_id <- raster::getValues(rc)    
  xy <- raster::xyFromCell(rc,1:raster::ncell(rc))
  df <- base::data.frame(xy, clump_id, is_clump = rc[] %in% raster::freq(rc, useNA = 'no')[,1])
  patch_coords <- plyr::ddply(df[df$is_clump == TRUE, ], plyr::.(clump_id), plyr::summarise, xm = base::mean(x), ym = base::mean(y))
  out <- base::list(patchrast=rc, coords=patch_coords)
  if(base::isTRUE(givedist)) {
    d <-  sp::spDists(sp::coordinates(patch_coords[,2:3]),longlat=TRUE)
    out <- base::c(out, base::list(distance=d))
  } 
}  

#' @rdname habitat
#' @export
is.habitat <- function(x) {
  inherits(x, "habitat")
}