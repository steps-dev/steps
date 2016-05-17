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
  if (isTRUE(readpoly)) 
    if (is.null(pypath)) {
    cmd <- Sys.which('OSGeo4W.bat')
    pypath <- 'gdal_polygonize'
    if(cmd=='') {
      cmd <- 'python'
      pypath <- Sys.which('gdal_polygonize.py')
      if (!file.exists(pypath)) 
        stop("Could not find gdal_polygonize.py or OSGeo4W on your system.") 
    }
  }
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists)) 
      stop(sprintf('File already exists: %s', 
                   toString(paste(outshape, c('shp', 'shx', 'dbf'), 
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    # require(raster)
    writeRaster(x, {f <- tempfile(fileext='.tif')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  
  system2(cmd, args=(
    sprintf('"%s" "%s" %s -f "ESRI Shapefile" "%s.shp"', 
            pypath, rastpath, ifelse(quietish, '-q ', ''), outshape)))
  
  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), 
                   verbose=!quietish)
    return(shp) 
  }
  return(NULL)
}