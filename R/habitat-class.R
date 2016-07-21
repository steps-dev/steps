#' @title habitat objects
#' @description Underlying habitat for dlmpr.
#' @rdname habitat
#' @name patches
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
#' @importFrom SDMTools PatchStat
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
#' foo <- patches(r2,distance=1000,p4s='+init=epsg:4283')

patches <- function(x, distance, p4s, givedist=TRUE) {
  if(!is(x, 'Raster')) x <- raster::raster(x)
  if(!is(p4s, 'CRS')) p4s <- sp::CRS(p4s)
  if(base::is.na(sp::proj4string(x))) stop(base::substitute(x), ' lacks a CRS.')
  mask <- x
  mask[mask == 0] <- NA
  ## create a buffer around original patch
  buff <- overlay(gridDistance(mask,origin=c(1),omit=NULL) <= distance, x, fun=function(a,b) ifelse(b > 0, b, a*2))
  rc <- raster::clump(buff,directions = 4)
  clump_id <- raster::getValues(rc)    
  xy <- raster::xyFromCell(rc,1:raster::ncell(rc))
  df <- base::data.frame(xy, clump_id, is_clump = rc[] %in% raster::freq(rc, useNA = 'no')[,1])
  patch_coords <- plyr::ddply(df[df$is_clump == TRUE, ], plyr::.(clump_id), plyr::summarise, xm = base::mean(x), ym = base::mean(y))
  out <- base::list(patchrast=rc, coords=patch_coords, patchstats = SDMTools::PatchStat(rc))
  if(base::isTRUE(givedist)) {
    d <-  sp::spDists(sp::coordinates(patch_coords[,2:3]),longlat=TRUE)
    out <- base::c(out, base::list(distance=d))
  } 
  base::class(out) <- "habitat"
  return(out)
}  

#' @rdname habitat
#' @export
is.patches <- function(x) {
  inherits(x, "patches")
}


## set up habitat with 
#' @rdname habitat
#' @name as.habitat
#' @param patches an object to turn into a \code{habitat} object. Currently
#'   this can either be a \link[raster]{raster}, a list or \code{NULL} (see \code{details}),
#'   though more approaches will be added in the future
#' @return an object of class \code{habitat}, essentially a dataframe
#'   containing the coordinates, area, population and features (as columns) for
#'   each patch (rows)
#' @author Nick Golding
#' @export
#' @examples
#' 
#' #'# create a habitat from a list containing a raster that represents a habitat suitability model / species distribution model.
#' habitat <- as.habitat(list(r2,population = as.population(t(rmultinom(1, 
#'                                size = 100, prob = c(0.8,0.2,0.01))))))
#' 
#' # create a default habitat
#' habitat <- as.habitat(NULL)
#'
#' # create a multipatch habitat
#' patches <- list(coordinates = data.frame(x=runif( 10, min=-10, max=10),
#'                                                     y=runif( 10, min=-10, max=10)),
#'                                area = as.data.frame(exp(-seq(.1,10,length.out = 10))*10),
#'                                population = as.population(t(rmultinom(1, 
#'                                size = 100, prob = c(0.8,0.2,0.01)))),
#'                                features = data.frame(temperature = 10))
#' habitat <- as.habitat(patches)
#'                                

as.habitat <- function (patches) {
  switch(class(patches)[1],
         NULL = habitatDefault(),
         list = list2habitat(patches))
}

#' @rdname habitat
#' @export
is.habitat <- function (x) inherits(x, 'habitat')

#' @rdname habitat
#' @param x an object to print or test as a habitat object
#' @param \dots further arguments passed to or from other methods.
#' @export
#' @examples
#' # print method
#' print(habitat)
#'
print.habitat <- function(x, ...) {
  if (is.null(nrow(x$distance))) text <- sprintf('an a-spatial habitat with %s patch\n',1)
  else text <- sprintf('a spatial habitat with %s patches\n',nrow(x$distance))
  cat(text)
}


#' @rdname habitat
#' @export
#' @param habitat an object of class \code{habitat}
#' @details the accessor functions \code{distance}, \code{area},
#'   \code{population} and \code{features} either return or set corresponding
#'   sub-dataframes of the \code{habitat} object
#' @examples
#' # get and set the area
#' area(habitat)
#'

area <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat$habitat[, attr(habitat$habitat, 'area'), drop = FALSE]
  ans <- squashhabitat(ans)
  return (ans)
}

#' @rdname habitat
#' @export
`area<-` <- function (habitat, value) {
  areaCheck(value)
  stopifnot(is.habitat(habitat))
  habitat$habitat[, attr(habitat$habitat, 'area')] <- value
  habitat
}

#' @rdname habitat
#' @name as.population
#' @param x a vector, data.frame or matrix of population(s) numbers for each stage(s) and each patch(s)
#' @export
#' @examples 
#' # starting population numbers for each step in the demographic function
#' population <- as.population(c(80,30,10))
#' population <- as.population(t(rmultinom(10, size = 100, prob = c(0.8,0.2,0.01))))

as.population <- function(x,...){
  if(base::is.null(base::dim(x))){
    names(x) <- base::paste0("stage",seq_along(x))
  } else { 
    x <- base::as.data.frame(x,...)
    ns <- base::dim(x)[[2]]
    np <- base::dim(x)[[1]]
    # s.names <- base::dimnames(x)[[2]]
    # p.names <- base::dimnames(x)[[1]]
    # if(base::is.null(s.names)) 
      s.names <- base::paste0("stage",1:ns)
    # if(base::is.null(p.names)) 
      p.names <- base::paste0("patch",1:np)
    base::dimnames(x) <- base::list(p.names, s.names)
  }
  base::class(x)<-c("population", class(x))
  return(x)
}

#' @rdname habitat
#' @name is.population
#' @export
#' @examples
#' population <- as.population(c(80,30,10))
#' is.population(population)
is.population <- function (x) {
  inherits(x, 'population')
}

#' @rdname habitat
#' @export
#' @examples
#'# get and set the population
#' population(habitat)
#' population(habitat) <- population(habitat) * 2
#' population(habitat)
#'
population <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat$habitat[, attr(habitat$habitat, 'population'), drop = FALSE]
  ans <- squashhabitat(ans)
  return (ans)
}

#' @rdname habitat
#' @export
`population<-` <- function (habitat, value) {
  stopifnot(is.habitat(habitat))
  populationCheck(value)
  stopifnot(all.equal(names(population(habitat)), names(value)))
  habitat$habitat[, attr(habitat$habitat, 'population')] <- value
  habitat
}


#' @rdname habitat
#' @name pop_patch_name
#' @param population matrix of states as cols and patches as rows.
pop_patch_name <- function (population) 
{
  states <- colnames(population)
  patches <- as.character(seq_len(nrow(population)))
  if (length(patches) == 1) {
    names <- states
  }
  else {
    names <- apply(expand.grid(states, patches), 1, paste, 
                   sep = "", collapse = "_patch_")
  }
  return(names)
}

#' @rdname habitat
#' @export
#' @examples
#'# get and set the features
#' features(habitat)

features <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat$habitat[, attr(habitat$habitat, 'features'), drop = FALSE]
  ans <- squashhabitat(ans)
  return (ans)
}

#' @rdname habitat
#' @export
`features<-` <- function (habitat, value) {
  stopifnot(is.habitat(habitat))
  stopifnot(is.data.frame(value))
  
  # for features, just overwrite whatever's there - including column numbers
  feature_cols <- attr(habitat$habitat, 'features')
  
  if (is.null(feature_cols) | length(feature_cols) == 0) {
    # if null (currently no features), add them
    attr(habitat$habitat, 'features') <- ncol(habitat) + 1:ncol(value)
  } else {
    # if not null (currently some features), overwrite them
    attrib <- attributes(habitat$habitat)
    attrib$names <- attrib$names[-feature_cols]
    # attrib$names
    habitat <- habitat$habitat[, -feature_cols]
    attributes(habitat) <- attrib
    attr(habitat, 'features') <- ncol(habitat) + seq_len(ncol(value))
  }
  habitat$habitat[, attr(habitat$habitat, 'features')] <- value
  habitat
}

#' @rdname habitat
#' @export
#' @examples
#'# get and set the features
#' spatial(habitat)
#'
spatial <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat$raster_patches
  return (ans)
}


#' @rdname habitat
#' @export
#' @examples
#'# get and set the features
#' coordinates(habitat)
#'
coordinates <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat$habitat[, attr(habitat$habitat, 'coordinates'), drop = FALSE]
  ans <- squashhabitat(ans)
  return (ans)
}

#' @rdname habitat
#' @export
`coordinates<-` <- function (habitat, value) {
  stopifnot(is.habitat(habitat))
  coordinatesCheck(value)
  stopifnot(all.equal(names(coordinates(habitat)), names(value)))
  habitat$habitat[, attr(habitat$habitat, 'coordinates')] <- value
  habitat
}

#' @rdname habitat
#' @export
#' @examples
#'# get and set the distance matrix
#' distance(habitat)
#'
distance <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat$distance
  return (ans)
}

#' @rdname habitat
#' @name is.distance
#' @export
is.distance <- function (x) {
  inherits(x, 'distance')
}

#' @rdname habitat
#' @param i index specifying the patches to include in the subset
#'   \code{habitat} object
#' @export
#' @examples
#' # habitats can be subsetted to get sub-habitats of patches with double
#' # braces
#' habitat
#' habitat[[1]]
#' habitat[[1:2]]
#'
`[[.habitat` <- function (x, i) {
  attrib <- attributes(x)
  attrib$row.names <- attrib$row.names[i]
  d <- attrib$distance[i, i, drop = FALSE]
  rownames(d) <- colnames(d) <- seq_along(i)
  attrib$distance <- d
  x <- squashhabitat(x)
  x <- x[i, ]
  attributes(x) <- attrib
  return (x)
}

#' 
#' @rdname habitat
#' @name K
#' @param x a raster of species habitit suitability (occupancy)
#' @param a parameter for carrying capacity model
#' @param b parameter for carrying capacity model
#' @param c parameter for carrying capacity model
#' @param type model form for converting occurrence to carrying capacity.
#' @author Skipton Woolley
K<-function(x,a=1,b=.2,c=0.5,type=c('exp','logit','linear')){
  type <- match.arg(type)
  switch(type,
         exp = exp((a*x)-b),
         linear = a*(x)-b,
         logit = a+(1/(1+exp(-b*x+c))))
}

raster2habitat <- function(input){ # will add in other options here later, but first lets get it working.
  # r[] <- scales::rescale(r[])
  r_id <- which(sapply(input,function(x)inherits(x,"RasterLayer")))
  r <- input[[r_id]]
  rthr <- r > stats::quantile(r[], .6) # default to .6 atm.
  patches <- dlmpr::patches(rthr, distance=1000, crs(rthr))  # default 1000m atm.
  
  if(!any(which(sapply(input,function(x)inherits(x,"population"))))){
  # lets calculated carrying capacity from occurrence.
  k <- K(r,type='exp') # will make this customisable one day.
  k_mask <- raster::mask(k,patches$patchrast)
  
  # calculate carrying capacity
  clump_id <- raster::getValues(patches$patchrast)    
  df <- base::data.frame(k=getValues(k_mask), clump_id, is_clump = patches$patchrast[] %in% raster::freq(patches$patchrast, useNA = 'no')[,1])
  intN <- plyr::ddply(df[df$is_clump == TRUE, ], plyr::.(clump_id), plyr::summarise, intN = base::sum(k))
  
  # calculate intial abundance per patch.
  intN <- intN[,2]*1.2 - (.1*patches$patchstats$perimeter)
  intN <- base::ifelse(intN<0,0,intN)
  
  # estimate population size based on starting K. # will need to replace this with stable states.
  population <- as.population(t(sapply(intN, function(x)rmultinom(1,size=x,prob=c(0.8,0.2,0.01))))) #replace these with stable states.
  } else {
  p_id <-  which(sapply(input,function(x)inherits(x,"population")))
  population <- as.data.frame(input[[p_id]])
  }
    # populations <- as.population(t(sapply(intN, function(x)rmultinom(1,size=x,prob=c(0.8,0.2,0.01)))))
  area <- patches$patchstats$area
  coordinates <- patches$coords[,-1]
  # populations <- data.frame(rep(population,nrow(coordinates)),nrow(coordinates),length(population))
  suppressWarnings(habitat <- list(habitat=   data.frame(coordinates,
                                              area = area,
                                              populations=population)))
  rownames(habitat$habitat) <- 1:nrow(habitat$habitat)
  
  # work out column numbers
  ncoord <- ncol(coordinates)
  narea <- 1
  npop <- ncol(population)
  
  attr(habitat$habitat, 'coordinates') <- seq_len(ncoord)
  attr(habitat$habitat, 'area') <- narea + ncoord
  attr(habitat$habitat, 'population') <- seq_len(npop) + narea + ncoord
  
  # add distance matrix
  habitat$distance <- as.matrix(dist(coordinates))
  attr(habitat$distance, 'distance') <- distance
  habitat$raster_patches <- patches$patchrast
  attr(habitat$raster_patches, 'spatial') <- raster_patches

  # set class
  class(habitat) <- c('habitat', class(habitat))
  
  # set class & return
  return (habitat)
}

list2habitat <- function (input) {
  
  # need to include the capacity to identify rasters (HSM)
  
  # check the elements
  if(any(sapply(input,function(x)inherits(x,"RasterLayer")))){
    habitat <- raster2habitat(input)
  } else {
         
  if(length(input)<4)stop()
    
  stopifnot(length(input) == 4)
  stopifnot(sort(names(input)) == c('area', 'coordinates', 'features', 'population'))
  stopifnot(all(sapply(input, is.data.frame)))
  
  # check components
  areaCheck(input$area)
  if(names(input$area)!='area')names(input$area)<-'area'
  populationCheck(input$population)
  
  # reset order and tidy up row names
  suppressWarnings(habitat <-  list(habitat=data.frame(input$coordinates,
                                           area = input$area,
                                           input$population,
                                           input$features)))
  rownames(habitat$habitat) <- 1:nrow(habitat$habitat)
  # work out column numbers
  ncoord <- ncol(input$coordinates)
  narea <- 1
  npop <- ncol(input$population)
  nfeat <- ncol(input$features)
  
  attr(habitat$habitat, 'coordinates') <- seq_len(ncoord)
  attr(habitat$habitat, 'area') <- narea + ncoord
  attr(habitat$habitat, 'population') <- seq_len(npop) + narea + ncoord
  attr(habitat$habitat, 'features') <- seq_len(nfeat) + npop + narea + ncoord
  
  # set class
  class(habitat) <- c('habitat', class(habitat))
  
  # add distance matrix
  habitat$distance <- as.matrix(dist(input$coordinates))
  }
  # set class & return
  return (habitat)
  
}

# default standalone habitat
habitatDefault <- function () {
  habitat_list <- list(coordinates = data.frame(x = 0, y = 0),
                         area = data.frame(area = 1),
                         population = data.frame()[1, ],
                         features = data.frame()[1, ])
  habitat <- list2habitat(habitat_list)
  return (habitat)
}

squashhabitat <- function (x) {
  # if an object is a habitat, remove the habitat class (to make it a
  # dataframe again)
  if (is.habitat(x)) {
    classes <- class(x)
    classes <- classes[-which(classes == 'habitat')]
    class(x) <- classes
  }
  return (x)
}

#internal function checks.
coordinateCheck <- function (coordinates) {2
  stopifnot(ncol(coordinates) == 2)
  stopifnot(is.numeric(coordinates[, 1]))
  stopifnot(all(is.finite(coordinates[, 1])))
}

areaCheck <- function (area) {
  stopifnot(ncol(area) == 1)
  stopifnot(is.numeric(area[, 1]))
  stopifnot(all(is.finite(area[, 1])))
  stopifnot(all(area[, 1] > 0))
}

populationCheck <- function (population) {
  stopifnot(all(sapply(population, is.finite)))
  stopifnot(all(sapply(population, function(x) all(x >= 0))))
}

distanceCheck <- function (distance, habitat) {
  stopifnot(is.matrix(distance))
  stopifnot(nrow(distance) == ncol(distance))
  stopifnot(nrow(distance) == nrow(habitat))
  stopifnot(all(is.finite(distance)))
  stopifnot(all(distance >= 0))
  stopifnot(all(diag(distance) == 0))
}
