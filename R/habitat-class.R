#' @title habitat objects
#' @description Habitat is an object that contains the spatial distribution of the populations, habitat suitability and carrying capacity for the landscape or seascape. Habitat requires either predefined rasters of population size for each life-history, habitat suitability map (e.g. a species distribution model) and carrying capacity. However, habitat suitability map is the only mandatory raster, population and carrying capacity can be provided as numeric values or functions which manipulate the habitat suitability map rasters to generate population per-cell and/or carrying capacity per-cell.
#' @rdname habitat
#' @name as.habitat
#' @param features A named list of landscapes (or seascape) features and parameters used for setting up the habitat for dynamic metapopulation models.
#' @details parameter details for habitat function.
#' \itemize{
#' \item{habitat_suitability_map}{ List a named list that contains a \link[raster]{raster} that represents habitat suitability for the species. This need to be probabilitic (between zero and one). Functions that manipulate the landscape will alter this throughout the dynamic metapopulation simulations. This can also be a \link[raster]{stack} or \link[raster]{brick}, if a raster stack or brick is provided, then then each raster in the stack or brick represents a temporal change in habitat suitability. For example, 10 rasters could represent habitat suitability over a ten year period, one per year.  
#' \item{populations}{ List starting populations for each life-history stage. Either a \link[sp]{SpatialPointsDataFrame} which has the population size for each life history linked to a coordinate within the extent of the habitat_suitability_map. A raster of population per-cell for each life-history stage, or finally a single integer of population size for each life-histroy stage.}
#' \item{carrying_capacity}{ List either a raster that represent the carrying capacity of adult populations for each cell; a function which manipulates the habitat_suitability_map and converts it to carrying capacity; or finally an integer which presents a maximum carrying capacity for all cells.}
#' }
#' @return an object of class \code{habitat}.
#' @author Nick Golding & Skip Woolley
#' @export
#' @examples
#' @examples
#' library(raster)
#' library(dhmpr)
#' set.seed(42)
#' xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
#' Dd <- as.matrix(dist(xy))
#' w <- exp(-1/nrow(xy) * Dd)
#' Ww <- chol(w)
#' xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
#' coordinates(xy) <- ~x+y
#' r <- rasterize(xy, raster(points2grid(xy)), 'z')
#' proj4string(r) <- '+init=epsg:4283'
#' r[] <- scales::rescale(r[],to=c(0,1))
#' 
#' ## create a habitat from a list containing a habitat suitability raster and numeric values for population and carrying capacity.
#' features <- list('habitat_suitability_map'=r,
#'                        'population'=c(80,20,10),
#'                        'carrying_capacity'=100)
#' habs <- as.habitat(features)
#'                        
#' ## create a habitat from a list containing a habitat suitability raster, a SpatialPointsDataFrame for population and numeric values carrying capacity.
#' random_populations <- sampleRandom(r, size=50, na.rm=TRUE, sp=TRUE) 
#' random_populations@data <- as.data.frame(t(rmultinom(50, size = 100, prob = c(0.8,0.2,0.1))))
#' features <- list('habitat_suitability_map'=r,
#'                        'population'=random_populations,
#'                        'carrying_capacity'=100)
#'                                                 
#' habs <- as.habitat(features)
#'                                
#' #create a habitat from a list containing just a species distribution model will estimate populations per-patch.
#' habs <- as.habitat(list(r))

as.habitat <- function (features,...) {
           stopifnot(is.list(features))  
           if(is.list(features))list2habitat(features,...)
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
  text <- sprintf('a spatial habitat with %s cells (patches if you prefer)\n',length(x$habitat_suitability_map[!is.na(x$habitat_suitability_map[])]))
  cat(text)
}

#' @rdname habitat
#' @name as.habitat_suitability
#' @param x a raster, raster stack or raster brick of habitat suitability for the species.
#' @export
#' @examples 
#' # Underlying habitat suitability map
#' population <- as.habitat_suitability(r)

as.habitat_suitability <- function(x,...){
  stopifnot(inherits(x,c("RasterLayer","RasterBrick","RasterStack")))
  base::class(x)<-c("habitat_suitability", class(x))
  return(x)
}

#' @rdname habitat
#' @name is.habitat_suitability
#' @export
#' @examples
#' hsm <- as.habitat_suitability(r)
#' is.habitat_suitability(hsm)
is.habitat_suitability <- function (x) {
  inherits(x, 'habitat_suitability')
}

#' @rdname habitat
#' @name as.populations
#' @param x Starting populations for each life-history stage. Either a \link[sp]{SpatialPointsDataFrame} which has the population size for each life history linked to a coordinate within the extent of the habitat_suitability_map. A raster of population per-cell for each life-history stage, or finally a single integer of population size for each life-histroy stage.
#' @export
#' @examples 
#' # Underlying habitat suitability map
#' random_populations <- sampleRandom(r, size=50, na.rm=TRUE, sp=TRUE) 
#' random_populations@data <- as.data.frame(t(rmultinom(50, size = 100, prob = c(0.8,0.2,0.1))))
#' population <- as.populations(random_populations)

as.populations <- function(x,...){
  stopifnot(inherits(x,c("RasterBrick","RasterStack","SpatialPointsDataFrame","numeric")))
  base::class(x)<-c("populations", class(x))
  return(x)
}

#' @rdname habitat
#' @name is.populations
#' @export
#' @examples
#' is.populations(pops)
is.populations <- function (x) {
  inherits(x, 'populations')
}

#' @rdname habitat
#' @name as.population
#' @param x a vector, data.frame, matrix or raster of population(s) numbers for each stage(s) and each patch(s)
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
    if(base::is.null(s.names))
      s.names <- base::paste0("stage",1:ns)
    if(base::is.null(p.names))
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
  ans <- habitat$habitat_suitability_map
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
#' @param habitat
pop_patch_name <- function (habitat){
  states <- gsub('populations.','',colnames(population(habitat)))
  patches <- as.character(seq_len(nrow(population(habitat))))
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
#' get and set the features
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
#' get and set the features
#' patches(habitat)
patches <- function (habitat, which_stages=NULL) {
  stopifnot(is.habitat(habitat))
  if(is.null(stage)) which_stages <- stages(habitat) 
  pops <- population(habitat)[[which_stages]]
  ans <- data.frame(patch_id=cellFromXY(pops,rasterToPoints(pops)[,1:2]),population=rasterToPoints(pops)[,-1:-2])
  return (ans)
}

#' @rdname habitat
#' @export
#' @examples
#' # get and set the features
#' habitat_suitability(habitat)

habitat_suitability <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ans <- habitat$habitat_suitability_map
  return (ans)
}

#' @rdname habitat
#' @export
`habitat_suitability<-` <- function (habitat, updated_habitat_suitability_map) {
  stopifnot(is.habitat(habitat))
  suitabilityCheck(updated_habitat_suitability_map)
  habitat$habitat_suitability_map <- updated_habitat_suitability_map
  habitat
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
  dist <- x$distance
  hab <- x$habitat
  attrib_habitat <- attributes(hab)
  attrib_habitat$row.names <- attrib_habitat$row.names[i]
  d <- dist[i, i, drop = FALSE]
  rownames(d) <- colnames(d) <- seq_along(i)
  x$habitat <- squashhabitat(x$habitat)
  x$habitat <- x$habitat[i, ]
  attributes(x$habitat) <- attrib_habitat
  return (x)
}

#' 
#' @rdname habitat
#' @name  carrying_capacity_function
#' @param x a raster of species habitit suitability (occupancy).
#' @param type model form for converting occurrence to carrying capacity.
#' @param list parameters used to convert a habitat suitability map to carrying capacity. 
#' @author Skipton Woolley
carrying_capacity_function <- function(x,type=c('exp','logit','linear','custom'),...){
  print(as.list(match.call(x)))
  type <- match.arg(type)
  switch(type,
         exp = exp((a*x)-b),
         linear = a*(x)-b,
         logit = a+(1/(1+exp(-b*x+c))),
         custom = custom_fun)
}


list2habitat <- function (input) {

  # need to include the capacity to identify rasters (HSM)
  # check the elements
  which_inputs_are_rasters <- which(sapply(input,function(x)inherits(x,"RasterLayer")))
  which_inputs_are_not_rasters <- which(sapply(input,function(x)!inherits(x,"RasterLayer")))




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

raster2habitat <- function(input){ # will add in other options here later, but first lets get it working.
  r_id <- which(sapply(input,function(x)inherits(x,"RasterLayer")))
  habitat_suitability <- input[[r_id]]

  if(!any(which(sapply(input,function(x)inherits(x,"population"))))){
  # lets calculated carrying capacity from habitat suitability.
  initial_k <- carrying_capacity(habitat_suitability, 
                                 type=input$carrying_capacity_model,
                                 input$carrying_capacity_params)
  
  # estimate population size based on starting K. # will need to replace this with stable states.
  population <- as.population(t(sapply(initial_k, function(x)rmultinom(1,size=x,prob=c(0.8,0.2,0.01))))) #replace these with stable states.
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
  habitat$starting_raster <- r
  attr(habitat$raster_patches, 'patches') <- habitat$raster_patches
  attr(habitat$starting_raster, 'habitat_suitability') <- habitat$starting_raster

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

patchCheck <- function(spatial_patches){
  stopifnot(class(spatial_patches)=="RasterLayer")
}

suitabilityCheck <- function(suitability){
  stopifnot(class(suitability)=="RasterLayer")
}

assign