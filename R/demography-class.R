#' @title demography objects 
#' @name as.demography
#' @rdname demography
#' @description demography is one of the main functions for dhmpr, it helps you construct and process stage based matricies.
#' Once a demography object is created \code{summary} will return:
#'  The finite rate of increase ("lambda"),
#'  the stable stage distribution,
#'  the reproductive value,
#'  and the sensitivities and elasticities matrices.
#'  \code{plot} will return a graph object that plots the transitions between and amoungest each stage of the population matrix.
#' @param x For as.demography, x is a square matrix, that has transition stages between population stages.
#' @param ... other function calls.
#' @return An object of class demography, i.e, resulting from as.demography.
#' @author Skipton Woolley
#' @export
#' @examples 
#' mat <- matrix(c(.53,0,.52,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult') 
#' trans <- as.demography(mat)

as.demography <- function(x, ...){
    object <- list(...)
    if(length(object)==0)object <- list(1)
    if(sapply(object, is.dispersal)) d <- TRUE
    else d <- FALSE
    x <- base::as.matrix(x)
    if(base::diff(base::dim(x)) !=0) stop("Needs to be a square matrix with transition probabilities between each stage.")
    stage_matrixCheck(x) 
    di <- base::dim(x)[[1]]
    m.names <- base::dimnames(x)[[1]]
    if(base::is.null(m.names)) m.names <- base::paste0("stage.",1:di)
    base::dimnames(x) <- base::list(m.names, m.names)
    if (d) {
      check_trans_disp(x,object)
      demography <- structure(list(stage_matrix=x,dispersal=object[[1]]),class='demography')
    }
    else demography <- structure(list(stage_matrix=x),class='demography')
    if (!is.demography(demography)) {
      class(demography) <- c('demography', class(demography))
    }
    return(demography)
}

#' @rdname demography
#' @param object an object of \code{demography} class
#' @export
#' @examples
#' summary(tmat)
#' 
summary.demography <-
  function(x,...){
    demographyCheck(x)
    x <- x$stage_matrix
    name.mat<-deparse(substitute(x))
    x <- x
    di <- base::dim(x)[1]
    m.names <- base::dimnames(x)[[1]] 
    ea<- base::eigen(x)
    lambda <- base::abs(ea$values[1]) 
    ssd <- base::abs(ea$vectors[,1]/base::sum(ea$vectors[,1]) ) 
    ae <- base::eigen(base::t(x))
    vr <- base::abs(ae$vectors[,1]/ae$vectors[1,1] )
    sensitivity <-  (vr  %*%  base::t(ssd))  / (base::t(vr) %*% ssd)[1,1]
    elasticity <- sensitivity * x / lambda
    
    result<- list(lambda=lambda, stable.stage.distribution = ssd,
                  reproductive.value =vr, sensitivity = sensitivity,
                  elasticity=elasticity,name.mat=name.mat,m.names= m.names)
    return (result)
  }

#' @rdname demography
#' @name stages
#' @param demography a demography object
#' @export
#' @examples
#' # get component stages
#' stages(trans)

stages <- function (demography) {
  # given a list of demographys, extract all of the mentioned stages
  stages <- colnames(demography$stage_matrix)
  return (stages)
}

#' @rdname demography
#' @importFrom igraph graph.adjacency
#' @export
#' @author Nick Golding
#' @examples 
#' plot(tmat)
#' 

plot.demography <- function (x, ...) {
  # plot a dynamic using igraph
  
  # extract the transition matrix & create an igraph graph x
  par(mar=c(2,4,4,2))
  demographyCheck(x)
  x <- x$stage_matrix
  textmat <- base::t(x)
  textmat[textmat>0]<-base::paste0('p(',base::as.character(textmat[textmat>0]),')')
  textmat[textmat=='0'] <-''
  linkmat <- textmat != ''
  g <- igraph::graph.adjacency(linkmat, weighted = TRUE)
  
  # extract edge labels
  labels <- textmat[igraph::get.edges(g, base::seq_len(base::sum(linkmat)))]
  
  # vertex plotting details
  igraph::V(g)$color <- grDevices::grey(0.9)
  igraph::V(g)$label.color <- grDevices::grey(0.4)
  igraph::V(g)$label.family <- 'sans'
  igraph::V(g)$size <- 50
  igraph::V(g)$frame.color <- NA
  
  # edge plotting details
  igraph::E(g)$color <- grDevices::grey(0.5)
  igraph::E(g)$curved <- igraph::curve_multiple(g, 0.1)
  igraph::E(g)$arrow.size <- 0.5
  igraph::E(g)$label <- labels
  igraph::E(g)$loop.angle <- 4
  igraph::E(g)$label.color <- grDevices::grey(0.4)
  
  graphics::plot(g)
  
  # return the igraph x
  return (base::invisible(g))
  
}

#' @rdname demography
#' @name is.demography
#' @export
#' @examples
#' is.demography(tmat)
is.demography <- function (x) {
  inherits(x, 'demography')
}


stage_matrixCheck <- function (x) {
  stopifnot(ncol(x) == nrow(x))
  stopifnot(is.matrix(x))
  stopifnot(all(is.finite(x)))
}

demographyCheck <- function (x) {
  stopifnot(inherits(x,'demography'))
  stopifnot(is.list(x))
}

dots <- function(...) {
  eval(substitute(alist(...)))
}


#### Casey I realised I messed up the demographic projections.
#### I had the matrix multipication around the wrong way.
#### I had vec_pops%*%stage_matrix, where is should have been stage_matrix%*%vec_pops.
#### I'm going to write a function which shouild do all the demographic projections 
#### and hopefully sort out any issues.

#'@param demographic_params a list of parameters which can be used to manipulate demographic projections.
#'\itemize{
#'#'\item{stage_matrix_sd}{Matrix with the standard deviation of the probabilities in \code{mat}.}
#'}

estimate_demography <- function(demography, populations, parameters){
   
  pop_vec <- lapply(populations(habitat),function(x)c(x[]))
  pop_mat <- do.call(cbind,pop_vec)
  ns <- length(stages(demography))
  
   if(all(dim(demography$stage_matrix)==ns)){
     # message for Casey: This will estimate deterministic population growth for a global transition matrix 
     
      message('Using a global demographic transition matrix for all patches\n')
      pop_mat_new <- t(sapply(1:nrow(pop_mat),function(x)estdemo(pop_mat[x,],demography$stage_matrix))) 
  } else {
    # message for Casey: This will estimate deterministic population growth for a local (one percell). 
    # The demographic$stage_matrix might need to be renamed to global_stage_matrix and local_stage_matricies. 
    # local_stage_matrices should have the dimensions ncells(rows),nstages*nstages(cols).
    
      message('Using a local demographic transition matricies; one per patch\n')
      pop_mat_new <- t(sapply(1:nrow(pop_mat),function(x)estdemo(pop_mat[x,],matrix(demography$stage_matrix[x,],ns,ns))))
  }
  # might need to return this as populations so we can slot back into: populations(habitat) <- new_populations.
  return(pop_mat_new)
}


estdemo <- function(popvec, stage_matrix, stage_matrix_sd=NULL){
  
    # no stochatisity
    popvec_new <- stage_matrix%*%popvec
    return(popvec_new)
  
}


