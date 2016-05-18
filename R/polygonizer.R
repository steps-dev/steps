#' polygonizer
#' Convert raster data to a ESRI polygon shapefile and (optionally) a SpatialPolygonsDataFrame
#' 
#' @param x an R Raster layer, or the file path to a raster file recognised by GDAL 
#' @param outshape the path to the output shapefile (if NULL, a temporary file will 
#'           be created) 
#' @param pypath the path to gdal_polygonize.py or OSGeo4W.bat (if NULL, the function 
#'        will attempt to determine the location)
#' @param readpoly should the polygon shapefile be read back into R, and returned by
#'           this function? (logical) 
#' @param quietish should (some) messages be suppressed? (logical)
#' @author John Baumgartner
#' @export

polygonizer <- function(x, outshape=NULL, pypath=NULL, readpoly=TRUE, 
                        quietish=TRUE) {
  if (base::isTRUE(readpoly)) 
    if (base::is.null(pypath)) {
    cmd <- base::Sys.which('OSGeo4W.bat')
    pypath <- 'gdal_polygonize'
    if(cmd=='') {
      cmd <- 'python'
      pypath <- base::Sys.which('gdal_polygonize.py')
      if (!base::file.exists(pypath)) 
        stop("Could not find gdal_polygonize.py or OSGeo4W on your system.") 
    }
  }
  if (!base::is.null(outshape)) {
    outshape <- base::sub('\\.shp$', '', outshape)
    f.exists <- base::file.exists(base::paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (base::any(f.exists)) 
      stop(base::sprintf('File already exists: %s', 
                   toString(base::paste(outshape, c('shp', 'shx', 'dbf'), 
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- base::tempfile()
  if (is(x, 'Raster')) {
    raster::writeRaster(x, {f <- base::tempfile(fileext='.tif')})
    rastpath <- base::normalizePath(f)
  } else if (base::is.character(x)) {
    rastpath <- base::normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  
  base::system2(cmd, args=(
    base::sprintf('"%s" "%s" %s -f "ESRI Shapefile" "%s.shp"', 
            pypath, rastpath, base::ifelse(quietish, '-q ', ''), outshape)))
  
  if (base::isTRUE(readpoly)) {
    shp <- rgdal::readOGR(base::dirname(outshape), layer = base::basename(outshape), 
                   verbose=!quietish)
    return(shp) 
  }
  return(NULL)
}