#'hsm_to_metapop
#'
#' @param hsm a occupancy model Raster layer
#' @param distance the neighbourhood distance. Patches that occur within this distance of
#'  one another will be clumped. This should be in the units of the CRS given
#'  in p4s. If this is zero, the function will identify patches defined by 
#'  8-connectivity (i.e. queen's case).
#' @param p4s an equal-area projection appropriate for the region. This will be used when 
#'      calculating inter-patch distances. The returned objects will be in the 
#'      original coordinate system.
#' @param givedist should the distance matrix be returned? (logical). Distances are in the 
#'           units of p4s, and are shortest cartesian distances between patches.
#' @param x a raster of species habitit suitability (occupancy)
#' @param a parameter for carrying capacity model
#' @param b parameter for carrying capacity model
#' @param c parameter for carrying capacity model
#' @param cc_mod model form for converting hs to carrying capacity. 

#' @importFrom raster raster
#' @export
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
#' metapop_spatial <- hsm_to_metapop(r2,occ_thresh=0.8, distance=1000, p4s='+init=epsg:3577', 
#'                                   givedist=TRUE, a=30,b=2,c=0.5,cc_mod='exp')
hsm_to_metapop <- function(hsm, occ_thresh=0.8, distance, p4s, 
                           givedist=TRUE, a=6,b=3,c=0.5,cc_mod='exp'){
  rthr <- hsm > stats::quantile(hsm[], occ_thresh)
  patches <- dlmpr::patchify(rthr,distance = distance,p4s=p4s)  
  K <- dlmpr::cc_fun(hsm,a=a,b=b,c=c,cc_mod=cc_mod)
  totalHS <- raster::extract(hsm,patches$patchpoly,sum,na.rm=TRUE)
  meanHS <- raster::extract(hsm,patches$patchpoly,mean,na.rm=TRUE)
  patchK <- raster::extract(K,patches$patchpoly,sum,na.rm=TRUE)
  ppgs <- patches$patchpoly
  pp <- methods::slot(ppgs, "data")
  di <- base::unlist(base::lapply(1:base::nrow(pp), function(i) {
    x <- sp::coordinates(ppgs@polygons[[i]]@Polygons[[1]])
    y <- base::rbind(x[-1,], x[1,])
    d <- raster::pointDistance(x, y, lonlat=FALSE)
    perimeter <- base::sum(d)
    perimeter
    }))
  intN <- raster::extract(K,patches$patchpoly,fun=function(x,...)base::sum(!base::is.na(x)))
  intN <- intN*1.2 - (.1*di)
  intN <- base::ifelse(intN<0,0,intN)
  mean_coords <- sp::getSpPPolygonsLabptSlots(ppgs)
  base::colnames(mean_coords) <- c('x','y')
  results <- base::list(hsm=hsm,patches=patches,k=K,
                  indices=base::data.frame(totalHS=totalHS,meanHS=meanHS,patchK=patchK,intN=intN),
                  mean_coords=mean_coords)
  return(results)
}