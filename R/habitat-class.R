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
#' hsm <- as.habitat_suitability(r)
#' pops <- as.populations(c(80,20,10))
#' cc <- as.carrying_capacity(100)
#' 
#' features <- list('habitat_suitability_map'=hsm,
#'                  'population'=pops,
#'                  'carrying_capacity'=cc)
#' habs <- as.habitat(features)
#'                        
#' ## create a habitat from a list containing a habitat suitability raster, a SpatialPointsDataFrame for population and numeric values carrying capacity.
#' random_populations <- sampleRandom(r, size=50, na.rm=TRUE, sp=TRUE) 
#' random_populations@data <- as.data.frame(t(rmultinom(50, size = 100, prob = c(0.8,0.2,0.1))))
#' features <- list('habitat_suitability_map'=as.habitat_suitability(r),
#'                        'population'=as.populations(random_populations),
#'                        'carrying_capacity'=as.carrying_capacity(100))
#'                                                 
#' habs <- as.habitat(features)

#########################
### habitat functions ###
#########################

as.habitat <- function (features,...) {
           stopifnot(is.list(features))
           stopifnot(length(features)==3)
           if(!is.habitat_suitability(features[[1]]))stop('first object in list must be "habitat_suitability"')
           if(!is.populations(features[[2]]))stop('second object in list must be "populations"')
           if(!is.carrying_capacity(features[[3]]))stop('third object in list must be "carrying_capacity"')
           features  
           # transformed_habitat <- list2habitat(features}
           # return(transformed_habitat)
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

#####################################
### habitat suitability functions ###
#####################################

#' @rdname habitat
#' @name as.habitat_suitability
#' @param x a raster, raster stack or raster brick of habitat suitability for the species.
#' @export
#' @examples 
#' # Underlying habitat suitability map
#' hsm <- as.habitat_suitability(r)

as.habitat_suitability <- function(x,...){
  stopifnot(inherits(x,c("RasterLayer","RasterBrick","RasterStack")))
  attr(x, "habitat") <- "habitat_suitability"
  # this causes conflicts with s4 classes.
  # base::class(x)<-c("habitat_suitability", class(x))
  return(x)
}

#' @rdname habitat
#' @name is.habitat_suitability
#' @export
#' @examples
#' hsm <- as.habitat_suitability(r)
#' is.habitat_suitability(hsm)
is.habitat_suitability <- function (x) {
  attr(x, 'habitat')=="habitat_suitability"
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

############################
### population functions ###
############################

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
  stopifnot(inherits(x,c("RasterLayer","RasterBrick","RasterStack","SpatialPointsDataFrame","numeric")))
  attr(x, "habitat") <- "populations"
  # base::class(x)<-c("populations", class(x))
  return(x)
}

#' @rdname habitat
#' @name is.populations
#' @export
#' @examples
#' is.populations(pops)
is.populations <- function (x) {
  attr(x, 'habitat')=="populations"
}

#' @rdname habitat
#' @export
#' @examples
#' # get and set the population
#' population(habitat)
#' population(habitat) <- population(habitat) * 2
#' population(habitat)

population <- function (habitat) {
  stopifnot(is.habitat(habitat))
  # habitat[[which(sapply(habitat,attr,"habitat")=='populations')]]
  pop <- habitat[[which(sapply(habitat,attr,"habitat")=='populations')]]
  # ans <- squashhabitat(ans)
  return (pop)
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
#' patches(habitat)
patches <- function (habitat, which_stages=NULL) {
  stopifnot(is.habitat(habitat))
  if(is.null(stage)) which_stages <- stages(habitat) 
  pops <- population(habitat)[[which_stages]]
  ans <- data.frame(patch_id=cellFromXY(pops,rasterToPoints(pops)[,1:2]),population=rasterToPoints(pops)[,-1:-2])
  return (ans)
}

###################################
### Carrying capacity functions ###
###################################

#' @rdname habitat
#' @name as.carrying_capacity
#' @export
#' @examples 


as.carrying_capacity <- function(x,...){
  stopifnot(inherits(x,c("RasterLayer","RasterBrick","RasterStack","SpatialPointsDataFrame","numeric","function")))
  attr(x, "habitat") <- "carrying_capacity"
  # base::class(x)<-c("populations", class(x))
  return(x)
}

#' @rdname habitat
#' @name is.carrying_capacity
#' @export
#' @examples
#' is.populations(pops)
is.carrying_capacity <- function (x) {
  attr(x, 'habitat')=="carrying_capacity"
}

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


##################################
### internal habitat functions ###
##################################

# this function will set up all the required needs of the experiment step.
# rasters will be left as rasters.
# if habitat suitability is a stack or brick, this infers that they represent a temporal change in habitat suitability.


list2habitat <- function (input) {

  # get the relevant habitat features: input[[1]] = habitat suitability.
  which_inputs_are_rasters <- which(sapply(input,function(x)inherits(x,c("RasterLayer"))))
  which_inputs_are_stacks <- which(sapply(input,function(x)inherits(x,c("RasterStack","RasterBrick"))))
  which_inputs_are_not_rasters <- which(sapply(input,function(x)!inherits(x,c("RasterLayer","RasterStack","RasterBrick"))))

  # generate an empty list to populate with habitat elements.
  habitat_list <- list()
  
  # sort out the habitat suitability maps.
  nrasters <- 0
  if(inherits(input[[1]],c("RasterStack","RasterBrick"))){
    nhsm <- nlayers(input[[1]])
    names(input[[1]]) <- paste0("habitat_suitability_t",1:nhsm)
    habitat_list <- sapply(1:nhsm,function(x){habitat_list[[x]]<-input[[1]][[x]]})
  } else {
    nhsm <- 1
    names(input[[1]]) <- "habitat_suitability"
    habitat_list[[1]] <- input[[1]]
  }
  # keep track of the number of rasters in  habitat.
  nrasters <- nrasters + nhsm
  
  # sort out the populations inputs.
  # if population is not a raster generate a raster from numeric or spatial points data frame.
  if(!inherits(inputs[[2]],c("RasterLayer","RasterStack","RasterBrick"))){
    inputs[[2]] <- populations_to_raster(input[[2]])
  }
  
  if(inherits(input[[2]],c("RasterStack","RasterBrick"))){
    npops <- nlayers(input[[2]])
    habitat_list <- sapply(1:npops,function(x){habitat_list[[x+nrasters]]<-input[[2]][[x]]})
  } else {
    nhsm <- 1
    habitat_list[[1]] <- input[[1]]
  }
  
  # generate an area raster from habitat suitability.
  area <- area_of_region(input[[1]]) # previous checks means this is a raster of some kind. 
  area_check(area)
  attr(area, "habitat") <- "area"
  
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
  # set class & return
  return (habitat)
}

#internal function checks.
coordinate_check <- function (coordinates) {
  stopifnot(ncol(coordinates) == 2)
  stopifnot(is.numeric(coordinates[, 1]))
  stopifnot(all(is.finite(coordinates[, 1])))
}

area_check <- function (area) {
  stopifnot(inherits(area,'RasterLayer'))
  stopifnot(is.numeric(area[!is.na(area[])]))
  stopifnot(all(area[!is.na(area[])] > 0))
}

population_check <- function (population) {
  stopifnot(all(sapply(population, is.finite)))
  stopifnot(all(sapply(population, function(x) all(x >= 0))))
}

## internal function for converting non-raster populations into rasters. 
population2raster <- function(pops,hab_suit){
  if(!inherits(pops,c("SpatialPointsDataFrame","numeric","function")))stop()
  
  # if numeric make a raster that is has the same values everywhere.
  if(inherits(pops,'numeric')){
    n_stages <- length(pops)
    pop_list <- list()
    for(i in 1:n_stages){
      tmp <- hab_suit
      tmp[!is.na(tmp[])]<-pops[i]
      pop_list[[i]] <- tmp
    }
    pop_stack <- stack(pop_list)
    names(pop_stack)<- paste0('stage',1:n_stages)
  }
  
  # if spatial points data frame generate a raster stack from known populations.
  if(inherits(pops,'SpatialPointsDataFrame')){
    
  
  }
  
  # if spatial points data frame generate a raster stack from known populations.
  if(inherits(pops,'function')){
    
      
  }  
  return(pop_rast)
}

## area functions. 
## internal functions for estimating area percell - useful for density dependence ect. 
area_at_site <- function(study_area, site_coords){
  if(raster::isLonLat(study_area)){
    #calculate area based on area function
    #convert kms to ms
    area_rast <- raster::area(study_area)*1000
    area_study <- raster::mask(area_rast,study_area)
    site_area <- raster::extract(area_study,site_coords,na.rm=TRUE)
  } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    cell_area <- raster::res(study_area)[1]*raster::res(study_area)[2]
    site_area <- rep(cell_area,nrow(site_coords))
  }
  return(site_area)
}

area_of_region <- function(study_area){
  if(raster::isLonLat(study_area)){
    #calculate area based on area function convert kms to ms
    area_rast <- area(study_area)*1000
    area_study <- mask(area_rast,study_area)

  } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    area_of_cell <- xres(study_area) * yres(study_area)
    study_area[!is.na(study_area[])] <- area_of_cell
    area_study <- study_area
  }
  return(area_study)
}
