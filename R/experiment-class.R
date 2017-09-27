#' @title experiment object
#' @rdname experiment-class
#' @name experiment
#' @description code{experiment} is an object which stores all the relevent population dynamics, population size(s), habitat feature, dispersal parameters and modules and sets up a spatially explicit dynamic metapopulation models for \link[dhmpr]{simulation}.
#' 
#' @param dots that can contain:
#'  times List The number of time steps of the least frequent discrete event. e.g. yearly population growth/management.
#'  dynamics a population dynamics such as stage-based transition matices, See \link[dhmpr]{as.transition}.
#'  population a population object, see \link[dhmpr]{as.population}
#'  habitat a habitat object, see \link[dhmpr]{as.habitat}.
#'  dispersal a dispersal object, see \link[dhmpr]{as.dispersal}. 
#'  module a module object, see \link[dhmpr]{as.module}, which is used to manipulate the habitat.
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
#' pop_hab_experiment_w_dispersal <- experiment(trans,habs,disp)
#' 
#' ## create a module for manipulating the landscape see \link[dhmpr]{as.module} for more details.
#' ##Create a named list with corresponding parameters and values
#' module_params = list(habitat=habs,
#'              fire_start_location = sample(ncell(suitability(habs)),10),
#'              prob = 0.24,
#'              continue_to_burn_prob = 0.01)
#'               
#' fire_module <- as.module(fire_spread,module_params)    
#' 
#' ## add module 
#' pop_hab_disp_experiments_w_module <- experiment(trans,habs,disp,fire_module)

experiment <- function(...){

  objects <- capture_dots(...)
  
  ## things that we need for the experiment to work
  
  #once all the objects are captured, break them up into their respective classes.
  trans_ids <-  which(sapply(objects,is.transition))
  disp_ids <-  which(sapply(objects,is.dispersal))
  if (length(trans_ids)>1)trans_ids <- trans_ids[which(trans_ids!=disp_ids)]
  pop_ids <-  which(sapply(objects,is.population))
  habitat_ids <-  which(sapply(objects,is.habitat))
  module_ids <-  which(sapply(objects,is.module))
  
  experiment <- list()
  if(length(trans_ids)==1)
    experiment$transition <- objects[[trans_ids]]
  if(length(pop_ids)==1 & length(habitat_ids)!=1){
    experiment$population <- as.population(objects[[pop_ids]])
  } else {
    experiment$population <- population(objects[[habitat_ids]])
  }  
  if(length(habitat_ids)==1){
    experiment$habitat <- objects[[habitat_ids]]
  } else{
    experiment$habitat <- as.habitat(NULL) #if not habitat supplied, supply a one patch habitat.
  }    
  if(length(disp_ids)==1)
    experiment$dispersal <- objects[[disp_ids]]
  if(length(module_ids)==1)
    experiment$module <- objects[[module_ids]]
    
  experiment <- as.experiment(experiment)
  
  return(experiment)
}

#' @rdname experiment-class
#' @export
is.experiment <- function (x) inherits(x, 'experiment')

as.experiment <- function (x) {
  if (!is.experiment(x)) {
    class(x) <- c('experiment', class(x))
  }
  return (x)
}

#' @rdname experiment-class
#' @param x an object to print or test as a habitat object
#' @param \dots further arguments passed to or from other methods.
#' @export
#' @examples
#' # print method
#' print(population_experiments)
#'
print.experiment <- function(x, ...) {
  text <- sprintf('A experiment object that contains:\n %s ',paste0(attr(x,'names'),collapse = ' \n '))
  cat(text)
}


#' @rdname experiment-class
#' @param experiment an object of class \code{experiment}
#' @export
habitat <- function (experiment) {
  stopifnot(is.experiment(experiment))
  value <- experiment[['habitat']]
  return (value)
}

#' @rdname experiment-class
#' @param value an object of class \code{habitat} (for
#'   \code{habitat(experiment) <- value}) or the value to assign to the
#'   \code{distance}, \code{area}, \code{population}, or \code{features}
#'   elements of a \code{habitat} object

#' @export
`habitat<-` <- function (experiment, value) {
  stopifnot(is.experiment(experiment))
  stopifnot(is.habitat(value))
  attr(experiment, 'habitat') <- value
  return (experiment)
}

#' @rdname experiment-class
#' @export
transition <- function (experiment) {
  stopifnot(is.experiment(experiment))
  value <- experiment[['transition']]
  return (value)
}

#' @rdname experiment-class
#' @export
`transition<-` <- function (experiment, value) {
  stopifnot(is.experiment(experiment))
  stopifnot(is.transition(value))
  attr(experiment, 'transition') <- value
  return (experiment)
}

#' #' @rdname experiment-class
#' #' @export
#' dispersal <- function (experiment) {
#'   stopifnot(is.experiment(experiment))
#'   stopifnot(any(sapply(experiment,is.dispersal)))
#'   value <- experiment[['dispersal']]
#'   return (value)
#' }
#' 
#' #' @rdname experiment-class
#' #' @export
#' `dispersal<-` <- function (experiment, value) {
#'   stopifnot(is.experiment(experiment))
#'   stopifnot(any(sapply(experiment,is.dispersal)))
#'   attr(experiment, 'dispersal') <- value
#'   return (experiment)
#' }

#' @rdname experiment-class
#' @export
as.matrix.experiment <- function(experiment,...){

  # if patches = 1
  hab <- habitat(experiment)
  p <- nrow(coords(hab))

  #get the stage demography matrix

  if(!any(sapply(experiment,is.dispersal))){
    
    experiment_mat <- stage_demographic_matrix(experiment)
    
  } else {

  stage_demog_mat <- stage_demographic_matrix(experiment)
  stage_disp_mat <- stage_disperal_matrix(experiment)
  P <- make.P(experiment)
  #experiment matrix
  experiment_mat <-  Matrix::t(P) %*% stage_disp_mat %*% P %*% stage_demog_mat
  }

  return(experiment_mat)
}

stage_demographic_matrix <- function(experiment){
  # this function sets up the stage demographic matrix for the vec-matrix metapopulation matrix.
  # it takes a experiment object which contains patch, and stage matrix.

  #get transition object
  trans <- transition(experiment)

  #get the stage matrix
  mat <- trans$stage_matrix

  #get the habitat from experiment object
  hab <- habitat(experiment)

  #get the number of patches
  p <- nrow(coords(hab))

  #set up the kronecker mat as a sparse matrix
  stage_demog_mat <- Matrix::kronecker(Matrix::Diagonal(p),mat)

  return(stage_demog_mat)
}


stage_disperal_matrix <- function(experiment){#disp_probs,states,M_sub){

  #get transition object
  trans <- transition(experiment)

  #get the habitat from experiment object
  hab <- habitat(experiment)

  # get the number of patches
  p <- nrow(coords(hab))

  # get the states
  states <- states(trans)

  # get the number of states
  s <- length(states)

  #get dispersal from habitat
  disp <- dispersal(experiment)
  
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

make.P <- function(experiment){

  # get patches and stages from experiment object
  hab <- habitat(experiment)
  p <- nrow(coords(hab)) #(or/cells)
  s <- length(states(transition(experiment))) #stages

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