#' @title habitat objects
#' @name patches
#' @rdname habitat
#' @description Underlying habitat for dlmpr.
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
#' foo <- patches(rthr, distance=1000, '+init=epsg:3577')
#' plot(foo$patchrast)
#' text(foo$coords[, 2:3], labels = foo$coords[, 1])

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
  out <- base::list(patchrast=rc, coords=patch_coords)
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
#'   this can either be a dynamic, a list or \code{NULL} (see \code{details}),
#'   though more approaches will be added in the future
#' @return an object of class \code{habitat}, essentially a dataframe
#'   containing the coordinates, area, population and features (as columns) for
#'   each patch (rows)
#' @export
#' @examples

#' # create a default habitat
#' habitat <- as.habitat(NULL)
#'
#' # create a marginally more interesting one-patch habitat
#' habitat <- as.habitat(list(coordinates = data.frame(x=runif( 10, min=-10, max=10),
#'                                                     y=runif( 10, min=-10, max=10)),
#'                                area = as.data.frame(exp(-seq(.1,10,length.out = 10))*10),
#'                                population = as.population(t(rmultinom(10, 
#'                                size = 100, prob = c(0.8,0.2,0.01)))),
#'                                features = data.frame(temperature = 10)))
as.habitat <- function (patches) {
  switch(class(patches)[1],
         NULL = habitatDefault(),
         dynamic = dynamichabitatDefault(patches),
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
  text <- sprintf('habitat with %s patches\n',
                  nrow(x))
  cat(text)
}

#' @rdname habitat
#' @name as.distance
#' @param x a vector, data.frame or matrix of distance(s) numbers for each stage(s) and each patch(s)
#' @export
#' @examples 
#' # starting distance numbers for each step in the demographic function
#' distance <- as.distance(c(80,30,10))
#' distance <- as.distance(t(rmultinom(10, size = 100, prob = c(0.8,0.2,0.01))))

as.distance <- function(x,...){
  if(base::is.null(base::dim(x))){
    names(x) <- base::paste0("stage",seq_along(x))
  } else { 
    x <- base::as.data.frame(x,...)
    ns <- base::dim(x)[[2]]
    np <- base::dim(x)[[1]]
    s.names <- base::dimnames(x)[[2]]
    p.names <- base::dimnames(x)[[1]]
    if(base::is.null(s.names)) s.names <- base::paste0("stage",1:ns)
    if(base::is.null(p.names)) p.names <- base::paste0("patch",1:np)
    base::dimnames(x) <- base::list(p.names, s.names)
  }
  base::class(x)<-c("distance", class(x))
  return(x)
}

#' @rdname habitat
#' @name is.distance
#' @export
#' @examples
#' distance <- as.distance(c(80,30,10))
#' is.distance(distance)
is.distance <- function (x) {
  inherits(x, 'distance')
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
    s.names <- base::dimnames(x)[[2]]
    p.names <- base::dimnames(x)[[1]]
    if(base::is.null(s.names)) s.names <- base::paste0("stage",1:ns)
    if(base::is.null(p.names)) p.names <- base::paste0("patch",1:np)
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
  ans <- habitat[, attr(habitat, 'area'), drop = FALSE]
  ans <- squashhabitat(ans)
  return (ans)
}

#' @rdname habitat
#' @export
#' @examples
#'# get and set the population
#' population(habitat)

population <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat[, attr(habitat, 'population'), drop = FALSE]
  ans <- squashhabitat(ans)
  return (ans)
}

#' @rdname habitat
#' @export
#' @examples
#'# get and set the features
#' features(habitat)
#'
features <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat[, attr(habitat, 'features'), drop = FALSE]
  ans <- squashhabitat(ans)
  return (ans)
}

#' @rdname habitat
#' @export
#' @examples
#'# get and set the distance matrix
#' distance(habitat)
#'
distance <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- attr(habitat, 'distance')
  return (ans)
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

coordinates <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat[, attr(habitat, 'coordinates'), drop = FALSE]
  ans <- squashhabitat(ans)
  return (ans)
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

list2habitat <- function (list) {
  
  # check the elements
  stopifnot(length(list) == 4)
  stopifnot(sort(names(list)) == c('area', 'coordinates', 'features', 'population'))
  stopifnot(all(sapply(list, is.data.frame)))
  
  # check components
  areaCheck(list$area)
  populationCheck(list$population)
  
  # reset order and tidy up row names
  suppressWarnings(habitat <- data.frame(list$coordinates,
                                           area = list$area,
                                           list$population,
                                           list$features))
  rownames(habitat) <- 1:nrow(habitat)
  
  # work out column numbers
  ncoord <- ncol(list$coordinates)
  narea <- 1
  npop <- ncol(list$population)
  nfeat <- ncol(list$features)
  
  attr(habitat, 'coordinates') <- seq_len(ncoord)
  attr(habitat, 'area') <- narea + ncoord
  attr(habitat, 'population') <- seq_len(npop) + narea + ncoord
  attr(habitat, 'features') <- seq_len(nfeat) + npop + narea + ncoord
  
  # set class
  class(habitat) <- c('habitat', class(habitat))
  
  # add distance matrix
  coord <- coordinates(habitat)
  distance <- as.matrix(dist(coord))
  distance(habitat) <- distance
  
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

# default habitat for a dynamic
dynamichabitatDefault <- function (dynamic) {
  population <- as.list(rep(0, length(states(dynamic))))
  names(population) <- states(dynamic)
  population <- as.data.frame(population)
  habitat_list <- list(coordinates = data.frame(x = 0, y = 0),
                         area = data.frame(area = 1),
                         population = population,
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
