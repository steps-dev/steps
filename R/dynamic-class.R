#' @title dynamic object
#' @rdname dynamic 
#' @name 
#' @description code{dynamic} is an object which stores all the relevent 
#' transitions, population(s), habitat, dispersal and module objects for use in \link[dhmpr]{simulation} function.
#' The main difference between \link[pop]{dynamic} and \code{dynamic} is that dhmpr requires the input of set objects to run.
#' @param dots that can contain:
#'  transition a transition object, See \link[dhmpr]{as.transition}.
#'  population a population object, see \link[dhmpr]{as.population}
#'  habitat a habitat object, see \link[dhmpr]{as.habitat}.
#'  dispersal a dispersal object, see \link[dhmpr]{as.dispersal}. 
#'  module a module object, see \link[dhmpr]{as.module}.
#'  custom a function to manipulate transition, population, habitat, or dispersal see 'custom function'. \code{as.customfun}.
#' @export
#' @examples
#'
#' ## Create transition matrix
#' mat <- matrix(c(.53,0,.52,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult')
#' trans <- as.transition(mat)
#'
#' ## Create population
#' pop <- as.population(data.frame('larvae'=80,'juvenile'=29,'adult'=5) )
#'
#' ## Create habitat
#' library(raster)
#' set.seed(42)
#' xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
#' Dd <- as.matrix(dist(xy))
#' w <- exp(-1/nrow(xy) * Dd)
#' Ww <- chol(w)
#' xy$z <- t(Ww) %*% rnorm(nrow(xy), 0, 0.1)
#' coordinates(xy) <- ~x+y
#' r <- rasterize(xy, raster(points2grid(xy)), 'z')
#' hab <- raster(r)
#' res(hab) <- 0.01
#' hab <- resample(r, hab)
#' proj4string(hab) <- '+init=epsg:4283'
#' habs <- as.habitat(list(hab,population = pop))

#' ## Create dispersal
#' dispersal_params <- list(alpha=list('larvae'=2,'juvenile'=0,'adult'=3),
#'                probability=list('larvae'=0.2,'juvenile'=0,'adult'=0.6))
#'
#' disp <- as.dispersal(dispersal_params)
#' 
#' ## create a module for manipulating the landscape see \link[dhmpr]{as.module} for more details.
#' fun <- fire_spread
#'
#' ##Create a named list with corresponding parameters and values
#' module_params = list(habitat=habs,
#'              fire_start_location = sample(ncell(suitability(habs)),10),
#'              prob = 0.24,
#'              continue_to_burn_prob = 0.01)
#'               
#' fire_module <- as.module(fun,module_params)    
#' 
#' ## simple transition matrix and population as a dynamic object
#' population_dynamics <- dynamic(trans,pop)
#' 
#' ## simple transition matrix, population and habitat as a dynamic object
#' population_habitat_dynamics <- dynamic(trans,pop,habs) 
#' 
#' ## because habitat contrains population information we can drop pop
#' population_habitat_dynamics <- dynamic(trans,habs)
#' 
#' ## add dispersal 
#' pop_hab_dynamics_w_dispersal <- dynamic(trans,habs,disp)
#' 
#' ## add module 
#' pop_hab_disp_dynamics_w_module <- dynamic(trans,habs,disp,fire_module)

dynamic <- function(...){

  objects <- capture_dots(...)
  
  #once all the objects are captured, break them up into their respective classes.
  trans_ids <-  which(sapply(objects,is.transition))
  disp_ids <-  which(sapply(objects,is.dispersal))
  if (length(trans_ids)>1)trans_ids <- trans_ids[which(trans_ids!=disp_ids)]
  pop_ids <-  which(sapply(objects,is.population))
  habitat_ids <-  which(sapply(objects,is.habitat))
  module_ids <-  which(sapply(objects,is.module))
  
  dynamic <- list()
  if(length(trans_ids)==1)
    dynamic$transition <- objects[[trans_ids]]
  if(length(pop_ids)==1 & length(habitat_ids)!=1){
    dynamic$population <- as.population(objects[[pop_ids]])
  } else {
    dynamic$population <- population(objects[[habitat_ids]])
  }  
  if(length(habitat_ids)==1){
    dynamic$habitat <- objects[[habitat_ids]]
  } else{
    dynamic$habitat <- as.habitat(NULL) #if not habitat supplied, supply a one patch habitat.
  }    
  if(length(disp_ids)==1)
    dynamic$dispersal <- objects[[disp_ids]]
  if(length(module_ids)==1)
    dynamic$module <- objects[[module_ids]]
    
  dynamic <- as.dynamic(dynamic)
  
  return(dynamic)
}

#' @rdname dynamic
#' @export
is.dynamic <- function (x) inherits(x, 'dynamic')

as.dynamic <- function (x) {
  if (!is.dynamic(x)) {
    class(x) <- c('dynamic', class(x))
  }
  return (x)
}

#' @rdname dynamic
#' @param x an object to print or test as a habitat object
#' @param \dots further arguments passed to or from other methods.
#' @export
#' @examples
#' # print method
#' print(population_dynamics)
#'
print.dynamic <- function(x, ...) {
  text <- sprintf('A dynamic object that contains:\n %s ',paste0(attr(x,'names'),collapse = ' \n '))
  cat(text)
}


#' @rdname dynamic
#' @param dynamic an object of class \code{dynamic}
#' @export
habitat <- function (dynamic) {
  stopifnot(is.dynamic(dynamic))
  value <- dynamic[['habitat']]
  return (value)
}

#' @rdname dynamic
#' @param value an object of class \code{habitat} (for
#'   \code{habitat(dynamic) <- value}) or the value to assign to the
#'   \code{distance}, \code{area}, \code{population}, or \code{features}
#'   elements of a \code{habitat} object

#' @export
`habitat<-` <- function (dynamic, value) {
  stopifnot(is.dynamic(dynamic))
  stopifnot(is.habitat(value))
  attr(dynamic, 'habitat') <- value
  return (dynamic)
}

#' @rdname dynamic
#' @export
transition <- function (dynamic) {
  stopifnot(is.dynamic(dynamic))
  value <- dynamic[['transition']]
  return (value)
}

#' @rdname dynamic
#' @export
`transition<-` <- function (dynamic, value) {
  stopifnot(is.dynamic(dynamic))
  stopifnot(is.transition(value))
  attr(dynamic, 'transition') <- value
  return (dynamic)
}

#' @rdname dynamic
#' @export
dispersal <- function (dynamic) {
  stopifnot(is.dynamic(dynamic))
  stopifnot(any(sapply(dynamic,is.dispersal)))
  value <- dynamic[['dispersal']]
  return (value)
}

#' @rdname dynamic
#' @export
`dispersal<-` <- function (dynamic, value) {
  stopifnot(is.dynamic(dynamic))
  stopifnot(any(sapply(dynamic,is.dispersal)))
  attr(dynamic, 'dispersal') <- value
  return (dynamic)
}

#' @rdname dynamic
#' @export
as.matrix.dynamic <- function(dynamic,...){

  # if patches = 1
  hab <- habitat(dynamic)
  p <- nrow(coords(hab))

  #get the stage demography matrix

  if(!any(sapply(dynamic,is.dispersal))){
    
    dynamic_mat <- stage_demographic_matrix(dynamic)
    
  } else {

  stage_demog_mat <- stage_demographic_matrix(dynamic)
  stage_disp_mat <- stage_disperal_matrix(dynamic)
  P <- make.P(dynamic)
  #dynamic matrix
  dynamic_mat <-  Matrix::t(P) %*% stage_disp_mat %*% P %*% stage_demog_mat
  }

  return(dynamic_mat)
}

stage_demographic_matrix <- function(dynamic){
  # this function sets up the stage demographic matrix for the vec-matrix metapopulation matrix.
  # it takes a dynamic object which contains patch, and stage matrix.

  #get transition object
  trans <- transition(dynamic)

  #get the stage matrix
  mat <- trans$stage_matrix

  #get the habitat from dynamic object
  hab <- habitat(dynamic)

  #get the number of patches
  p <- nrow(coords(hab))

  #set up the kronecker mat as a sparse matrix
  stage_demog_mat <- Matrix::kronecker(Matrix::Diagonal(p),mat)

  return(stage_demog_mat)
}


stage_disperal_matrix <- function(dynamic){#disp_probs,states,M_sub){

  #get transition object
  trans <- transition(dynamic)

  #get the habitat from dynamic object
  hab <- habitat(dynamic)

  # get the number of patches
  p <- nrow(coords(hab))

  # get the states
  states <- states(trans)

  # get the number of states
  s <- length(states)

  #get dispersal from habitat
  disp <- dispersal(dynamic)
  
  # get dispersal
  disp_probs <- probdisp(disp,hab)

  # get the states that have dispersal matricies
  disp_states <- rep(0,s)
  disp_states[match(names(disp_probs),states)] <- match(names(disp_probs),states)

  # set up the sparse matrix
  dim_mat <-s*s*p^2
  stage_disp_mat <- Matrix::sparseMatrix(dims = c(1,dim_mat), i={}, j={},x=rep(0.0,dim_mat))

  #get the M-sub mat as a vec to save memory.
  dims<-p^2
  index<-0:(p-1)
  Index<-1+index*p+index
  M_sub<-sparseMatrix(i = rep(1,p), j = Index, x= rep(1.0,p), dims=c(1,dims))

  #Calculate the indices which have to be added for all p^2 indices of M_ii
  indx1 <- seq(1, p * p);
  xp <- (indx1 - 1) - ((indx1 -1) %/% p) * p
  yp <- (indx1 - 1) %/% p;
  indx2 <- xp + yp * s * p

  # add data to dispersal matrix
  z <- 1
  m <- 1
  for (i in seq_along(states)){
    if(i==disp_states[i]){   #check which stages have dispersal
      stage_disp_mat[z+indx2] <- disp_probs[[m]][indx1]
      m <- m + 1
    } else {
      stage_disp_mat[z+indx2] <- M_sub[indx1] # if no dispersal make it the M_sub mat.
    }
    z <- 1+(1+i-1)*(s*p^2+p)
  }
  dim(stage_disp_mat)<-rep(s*p,2)

  #check the resulting matrix
  # image(stage_disp_mat)

  return(stage_disp_mat)

}

make.P <- function(dynamic){

  # get patches and stages from dynamic object
  hab <- habitat(dynamic)
  p <- nrow(coords(hab)) #(or/cells)
  s <- length(states(transition(dynamic))) #stages

  rw<-cl<-rep(0, s*p)

  index<-0
  for( i in 1:s){
    for( j in 1:p){
      index<-index+1
      cl[index]<- (j-1)* s + i
      rw[index]<- (i-1)* p + j
    }
  }

  P <- Matrix::sparseMatrix(i = rw, j = cl, x= rep(1,length(rw)),dims=rep(s*p,2))
  return(P)
}

capture_dots <- function (...) {
  # capture arguments passed as dots, and grab their names even if not directly
  # named
  ans <- list(...)
  dots <- substitute(list(...))[-1]
  names(ans) <- sapply(dots, deparse)
  ans
}