#' @title habitat objects
#' @description Habitat is an object that contains the spatial distribution of the populations, habitat suitability and carrying capacity for the landscape or seascape. Habitat requires either predefined rasters of population size for each life-history, habitat suitability map (e.g. a species distribution model) and carrying capacity. However, habitat suitability map is the only mandatory raster, population and carrying capacity can be provided as numeric values or functions which manipulate the habitat suitability map rasters to generate population per-cell and/or carrying capacity per-cell.
#' @rdname habitat
#' @param features A named list of landscapes (or seascape) features and parameters used for setting up the habitat for dynamic meta-population models.
#' @details parameter details for habitat function.
#' \itemize{
#' \item{habitat_suitability_map}{ An object of class 'habitat_suitability', which must contain at least a single \link[raster]{raster} that represents habitat suitability for the species. This need to be probabilistic (between zero and one). Functions that manipulate the landscape will alter this throughout the dynamic meta-population simulations. This can also be a \link[raster]{stack} or \link[raster]{brick}, if a raster stack or brick is provided, then then each raster in the stack/brick represents habitat suitability at temporal time steps. For example, 10 rasters could represent habitat suitability over a ten year period, one per year.}  
#' \item{populations}{ List starting populations for each life-history stage. Either a \link[sp]{SpatialPointsDataFrame} which has the population size for each life history linked to a coordinate within the extent of the habitat_suitability_map. A raster of population per-cell for each life-history stage, or finally a single integer of population size for each life-histroy stage.}
#' \item{carrying_capacity}{ List either a raster that represent the carrying capacity of adult populations for each cell; a function which manipulates the habitat_suitability_map and converts it to carrying capacity; or finally an integer which presents a maximum carrying capacity for all cells.}
#' }
#' @return an object of class \code{habitat}.
#' @author Nick Golding & Skip Woolley
#' @export
#' @examples
#' 
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
#' cc <- as.carrying_capacity(300)
#' 
#' features <- list(hsm,pops,cc)
#' habitat <- as.habitat(features)
#'                        
#' ## create a habitat from a list containing a habitat suitability raster, a SpatialPointsDataFrame for population and numeric values carrying capacity.
#' random_populations <- sampleRandom(r, size=50, na.rm=TRUE, sp=TRUE) 
#' random_populations@data <- as.data.frame(t(rmultinom(50, size = 100, prob = c(0.8,0.2,0.1))))
#' features <- list('habitat_suitability_map'=as.habitat_suitability(r),
#'                  'population'=as.populations(random_populations),
#'                  'carrying_capacity'=as.carrying_capacity(100))
#'                                                 
#' habitat <- as.habitat(features)

#########################
### habitat functions ###
#########################

as.habitat <- function (features,...) {
           stopifnot(is.list(features))
           stopifnot(length(features)==3)
           if(!is.habitat_suitability(features[[1]]))stop('first object in list must be a "habitat_suitability" object')
           if(!is.populations(features[[2]]))stop('second object in list must be a "populations" object')
           if(!is.carrying_capacity(features[[3]]))stop('third object in list must be a "carrying_capacity" object')
           transformed_habitat <- list2habitat(features)
           return(transformed_habitat)
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

print.habitat <- function(x, ...) {
  numastext <- c('one','two','three','four','five','six','seven','eight','nine','ten')
  nhsms <- sum(lapply(x,attr,"habitat")=='habitat_suitability')
  nstages <- sum(lapply(x,attr,"habitat")=='populations')
  ncells <- length(x[[1]][!is.na(x[[1]])])
  if (all(c(nhsms,nstages)<11))  text <- sprintf('%s habitat suitability layer(s) with %s cells and \n%s population life-history stage(s).\n',numastext[nhsms],ncells,numastext[nstages])
  else text <- sprintf('%s habitat suitability layer(s) with %s cells and \n%s population life-history stage(s).\n',nhsms,ncells,nstages)
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
#'  
#' # Underlying habitat suitability map
#' hsm <- as.habitat_suitability(r)

as.habitat_suitability <- function(x,...){
  stopifnot(inherits(x,c("RasterLayer","RasterBrick","RasterStack")))
  attr(x, "habitat") <- "habitat_suitability"
  return(x)
}

#' @rdname habitat
#' @name is.habitat_suitability
#' @export
#' @examples
#' 
#' is.habitat_suitability(hsm)

is.habitat_suitability <- function (x) {
  attr(x, 'habitat')=="habitat_suitability"
}

#' @rdname habitat
#' @export
#' @examples
#' 
#' habitat_suitability(habitat)

habitat_suitability <- function (habitat,time_step=1) {
  stopifnot(is.habitat(habitat))
  ans <- habitat[which(sapply(habitat,attr,"habitat")=='habitat_suitability')][[time_step]]
  return (ans)
}

#' @rdname habitat
#' @export
`habitat_suitability<-` <- function (habitat, time_step=1, value) {
  stopifnot(is.habitat(habitat))
  # check_habitat_suitability(updated_habitat_suitability_map)
  habitat[which(sapply(habitat,attr,"habitat")=='habitat_suitability')][[time_step]] <- value
  habitat
}

############################
### population functions ###
############################

#' @rdname habitat
#' @name as.populations
#' @export
#' @examples
#' 
#' as.populations(c(80,60,20))

as.populations <- function(x,...){
  stopifnot(inherits(x,c("RasterLayer","RasterBrick","RasterStack","SpatialPointsDataFrame","numeric")))
  attr(x, "habitat") <- "populations"
  return(x)
}

#' @rdname habitat
#' @name is.populations
#' @export
#' @examples
#' 
#' is.populations(pops)
is.populations <- function (x) {
  attr(x, 'habitat')=="populations"
}

#' @rdname habitat
#' @export
#' @examples
#' 
#' # get and set the population
#' populations(habitat)

populations <- function (habitat, which_stages=NULL) {
  stopifnot(is.habitat(habitat))
  if(is.null(which_stages)) which_stages <- seq_len(sum(lapply(habitat,attr,"habitat")=='populations'))
  # habitat[[which(sapply(habitat,attr,"habitat")=='populations')]]
  pop <- habitat[which(sapply(habitat,attr,"habitat")=='populations')][which_stages]
  # ans <- squashhabitat(ans)
  return (pop)
}

#' @rdname habitat
#' @export
`populations<-` <- function (habitat, which_stages=NULL, value) { #which_stages=NULL,
  stopifnot(is.habitat(habitat))
  if(is.null(which_stages)) which_stages <- seq_len(sum(lapply(habitat,attr,"habitat")=='populations'))
  habitat[which(sapply(habitat,attr,"habitat")=='populations')][which_stages] <- value
  pops <- habitat
  return(pops)
}

###################################
### Carrying capacity functions ###
###################################

#' @rdname habitat
#' @name as.carrying_capacity
#' @export
#' @examples
#' as.carrying_capacity(100)

as.carrying_capacity <- function(x,...){
  stopifnot(inherits(x,c("RasterLayer","function")))
  attr(x, "habitat") <- "carrying_capacity"
  # base::class(x)<-c("populations", class(x))
  return(x)
}

#' @rdname habitat
#' @name is.carrying_capacity
#' @export
#' @examples
#' is.carrying_capacity(cc)
is.carrying_capacity <- function (x) {
  attr(x, 'habitat')=="carrying_capacity"
}

#' #' @rdname habitat
#' #' @name  carrying_capacity_function
#' #' @param x a raster of species habitat suitability (occupancy).
#' #' @param type model form for converting occurrence to carrying capacity.
#' #' @param list parameters used to convert a habitat suitability map to carrying capacity. 
#' 
#' carrying_capacity_function <- function(x, type=c('exp','logit','linear','custom'), custom_fun=NULL, params, ...){
#'   print(as.list(match.call(x)))
#'   type <- match.arg(type)
#'   switch(type,
#'     exp = exp((a*x)-b),
#'     linear = a*(x)-b,
#'     logit = a+(1/(1+exp(-b*x+c))),
#'     custom = custom_fun)
#' }

#' @rdname habitat
#' @export

carrying_capacity <- function (habitat,time_step=1) {
  stopifnot(is.habitat(habitat))
  # habitat[[which(sapply(habitat,attr,"habitat")=='carrying_capacity')]]
  cc <- habitat[which(sapply(habitat,attr,"habitat")=='carrying_capacity')][[time_step]]
  return(cc)
}

#' @rdname habitat
#' @export
`carrying_capacity<-` <- function (habitat, time_step=1, value) {
  stopifnot(is.habitat(habitat))
  # population_check(value)
  habitat[which(sapply(habitat,attr,"habitat")=='carrying_capacity')][[time_step]]<-value
  habitat
}

######################
### area functions ###
######################
#' @rdname habitat
#' @export

area <- function (habitat) {
  stopifnot(is.habitat(habitat))
  ar <- habitat[which(sapply(habitat,attr,"habitat")=='area')]
  return(ar)
}

## probably don't need this function as area of cells is unlikely to change.
#' #' @rdname habitat
#' #' @export
#' `area<-` <- function (habitat, new_area) {
#'   stopifnot(is.habitat(habitat))
#'   area_check(new_area)
#'   habitat[which(sapply(habitat,attr,"habitat")=='area')]<-new_area
#'   habitat
#' }

##################################
### internal habitat functions ###
##################################

# this function will set up all the required needs of the experiment step.
# rasters will be left as rasters.
# if habitat suitability is a stack or brick, this infers that they represent a temporal change in habitat suitability

list2habitat <- function (input) {

  # generate an empty list to populate with habitat elements.
  habitat_list <- list()
  
  # sort out the habitat suitability maps.
  nrasters <- 0
  if(inherits(input[[1]],c("RasterStack","RasterBrick"))){
    nhsm <- raster::nlayers(input[[1]])
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
  # update and make a brick.
  if(!inherits(input[[2]],c("RasterLayer","RasterStack","RasterBrick"))){
    input[[2]] <- populations2rasterbrick(input[[2]],input[[1]])
  }
  
  if(inherits(input[[2]],c("RasterStack","RasterBrick"))){
    npops <- nlayers(input[[2]])
    habitat_list[c(nrasters+1):c(nrasters+npops)] <- sapply(1:npops,function(x){habitat_list[[x+nrasters]]<-input[[2]][[x]]})
  } else {
    stop('check your population inputs, something is going wrong.')
  }
  
  # add npops to the total count of rasters. 
  nrasters <- nrasters + npops
  
  if(!inherits(input[[3]],c("RasterLayer","RasterBrick","RasterStack"))){
    input[[3]] <- carryingcapacity2raster(input[[3]],input[[1]])
  }
  
  # update the number of rasters
  nrasters <- nrasters + 1
  
  # add carrying capacity to the habitat list.
  habitat_list[[nrasters]] <- input[[3]]

  # generate an area raster from habitat suitability.
  nrasters <- nrasters + 1
  areas <- area_of_region(input[[1]]) # previous checks means this is a raster of some kind. 
  area_check(areas)
  habitat_list[[nrasters]] <- areas
  
  # This should correctly setup attributes. yay.
  habitat_list[1:nhsm] <- lapply(habitat_list[1:nhsm], `attr<-`, "habitat", "habitat_suitability")
  habitat_list[c(1+nhsm):c(nhsm+npops)] <- lapply(habitat_list[c(1+nhsm):c(nhsm+npops)], `attr<-`, "habitat", "populations")
  habitat_list[c(1+nhsm+npops)] <- lapply(habitat_list[c(1+nhsm+npops)], `attr<-`, "habitat", "carrying_capacity")
  habitat_list[c(2+nhsm+npops)] <- lapply(habitat_list[c(2+nhsm+npops)], `attr<-`, "habitat", "area")
  
  # set class
  class(habitat_list) <- c('habitat')

  # set class & return
  return (habitat_list)
}

#internal function checks.
area_check <- function (area) {
  stopifnot(inherits(area,c("RasterLayer","RasterBrick","RasterStack")))
  if(inherits(area,"RasterLayer")) {
    stopifnot(is.numeric(area[!is.na(area[])]))
    stopifnot(all(area[!is.na(area[])] > 0))    
  }else{
    stopifnot(sapply(area, function(x) is.numeric(x[!is.na(x[])])))
    stopifnot(sapply(area, function(x) all(x[!is.na(x[])] > 0)))  
  }
}

# population_check <- function (population) {
#   stopifnot(all(sapply(population, is.finite)))
#   stopifnot(all(sapply(population, function(x) all(x >= 0))))
# }

## internal function for converting non-raster populations into rasters. 
## raster bricks are more efficent. 
populations2rasterbrick <- function(pops,hab_suit,ss_dist){
  
  # I will in the capacity to alter population with a function soon.
  stopifnot(inherits(pops,c("SpatialPointsDataFrame","numeric"))) #,"function")))stop()
  if(inherits(hab_suit,c("RasterBrick","RasterStack"))) hab_suit <- hab_suit[[1]]
  
  # if numeric make a raster that is has the same values everywhere.
  if(inherits(pops,'numeric')){
    # if(length(pops) == 1 & ss_dist != NULL){
    #   cat('Initial populations were entered as a single number - they will be distributed \nacross the landscape based on the habitat suitability and stable age structure')
    #   pop_list <- list()
    #   if(length(ss_dist) != n_stages) stop('A stable stage muliplier needs to be entered for each life stage; please check your values.')
    #   for(i in 1:n_stages){
    #     tmp <- hab_suit
    #     tmp[!is.na(tmp[])] <- pops * ss_dist[i]
    #     pop_list[[i]] <- tmp * hab_suit
    #     }
    #   pop_brick <- raster::brick(pop_list)
    #   names(pop_brick)<- paste0('stage',1:n_stages)      
    # }else{
      n_stages <- length(pops)
      pop_list <- list()
      for(i in 1:n_stages){
        tmp <- hab_suit
        tmp[!is.na(tmp[])] <- pops[i]
        pop_list[[i]] <- tmp
      }
      pop_brick <- raster::brick(pop_list)
      names(pop_brick)<- paste0('stage',1:n_stages)     
    # }
  }
  
  # if spatial points data frame generate a raster stack from known populations.
  if(inherits(pops,'SpatialPointsDataFrame')){
    if(raster::projection(pops)!=raster::projection(hab_suit))stop('make sure your spatial points dataframe matches raster projections')
    n_stages <- ncol(pops@data)
    pop_brick <- raster::rasterize(pops,hab_suit)[[-1]]
    names(pop_brick) <- paste0('stage',1:n_stages)
    pop_brick <- raster::calc(pop_brick, fun=function(x){ x[is.na(x)] <- 0; return(x)})
    pop_brick <- raster::mask(pop_brick,hab_suit)
  }
 
  return(pop_brick)
}

## internal function for converting non-raster populations into rasters. 
## raster bricks are more efficent. 
carryingcapacity2raster <- function(carry_cap,hab_suit,...){
  
  # I will in the capacity to alter population with a function soon.
  stopifnot(inherits(carry_cap,c("numeric"))) #,"function")))stop()
  if(inherits(hab_suit,c("RasterBrick","RasterStack"))) hab_suit <- hab_suit[[1]]
  
  # if numeric make a raster that is has the same values everywhere.
  if(inherits(carry_cap,'numeric')){
      ccr <- hab_suit
      ccr[!is.na(ccr[])] <- carry_cap
      names(ccr)<- 'carrying_capacity'
  }
  return(ccr)
}


## area functions. 
## internal functions for estimating area percell - useful for density dependence etc. 
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
    area_rast <- raster::area(study_area)*1000
    area_study <- raster::mask(area_rast,study_area)

  } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    area_of_cell <- raster::xres(study_area) * raster::yres(study_area)
    study_area[!is.na(study_area[])] <- area_of_cell
    area_study <- study_area
  }
  return(area_study)
}
