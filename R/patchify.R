#' patchify
#' Clump raster cells into patches, with optional neighbourhood distance.
 
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
#' @importClassesFrom raster Raster
#' @importClassesFrom sp CRS

#' @examples
#' library(rasterVis)
#' library(raster)
#' set.seed(1)
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
#' levelplot(rthr, margin=FALSE, col.regions=0:1, colorkey=FALSE)
#' foo <- patchify(rthr, distance=1000, '+init=epsg:3577')
#' spplot(foo$patchrast, col.regions=rainbow(length(foo$patchpoly)))

patchify <- function(x, distance, p4s, givedist=TRUE) {
  if(!is(x, 'Raster')) x <- raster::raster(x)
  if(!is(p4s, 'CRS')) p4s <- sp::CRS(p4s)
  if(base::is.na(sp::proj4string(x))) stop(base::substitute(x), ' lacks a CRS.')
  x[x == 0] <- NA
  cc <- SDMTools::ConnCompLabel(x)
  p <- base::suppressWarnings(dlmpr::polygonizer(cc))
  sp::proj4string(p) <- sp::proj4string(x)
  pproj <- sp::spTransform(p, p4s)
  bproj <- rgeos::gBuffer(pproj, width=distance)
  pproj$patch <- sp::over(pproj, sp::disaggregate(bproj))
  b <- sp::spTransform(pproj, sp::CRS(sp::proj4string(x)))
  pout <- raster::aggregate(b, vars='patch')
  pout$patch <- base::as.factor(pout$patch)
  rout <- raster::rasterize(pout, x)
  out <- base::list(patchrast=rout, patchpoly=pout)
  if(base::isTRUE(givedist)) {
    poutproj <- sp::spTransform(pout, p4s)
    d <- rgeos::gDistance(poutproj, poutproj, byid=TRUE)
    out <- base::c(out, base::list(distance=d))
  } 
}  
