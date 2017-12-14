#' @title transition objects 
#' @name as.transition
#' @rdname transition
#' @description transition is one of the main functions for dhmpr, it helps you construct and process stage based matricies.
#' Once a transition object is created \code{summary} will return:
#'  The finite rate of increase ("lambda"),
#'  the stable stage distribution,
#'  the reproductive value,
#'  and the sensitivities and elasticities matrices.
#'  \code{plot} will return a graph object that plots the transitions between and amoungest each stage of the population matrix.
#' @param x For as.transition, x is a square matrix, that has transition states between population stages.
#' @param ... other function calls.
#' @return An object of class transition, i.e, resulting from as.transition.
#' @author Skipton Woolley
#' @export
#' @examples 
#' mat <- matrix(c(.53,0,.52,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' colnames(mat) <- rownames(mat) <- c('larvae','juvenile','adult') 
#' trans <- as.transition(mat)

as.transition <- function(x, ...){
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
      transition <- structure(list(stage_matrix=x,dispersal=object[[1]]),class='transition')
    }
    else transition <- structure(list(stage_matrix=x),class='transition')
    if (!is.transition(transition)) {
      class(transition) <- c('transition', class(transition))
    }
    return(transition)
}

#' @rdname transition
#' @param object an object of \code{transition} class
#' @export
#' @examples
#' summary(tmat)
#' 
summary.transition <-
  function(x,...){
    transitionCheck(x)
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

#' @rdname transition
#' @name states
#' @param transition a transition object
#' @export
#' @examples
#' # get component states
#' states(trans)

states <- function (transition) {
  # given a list of transitions, extract all of the mentioned states
  states <- colnames(transition$stage_matrix)
  return (states)
}

#' @rdname transition
#' @importFrom igraph graph.adjacency
#' @export
#' @author Nick Golding
#' @examples 
#' plot(tmat)
#' 

plot.transition <- function (x, ...) {
  # plot a dynamic using igraph
  
  # extract the transition matrix & create an igraph graph x
  par(mar=c(2,4,4,2))
  transitionCheck(x)
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

#' @rdname transition
#' @name is.transition
#' @export
#' @examples
#' is.transition(tmat)
#' 
is.transition <- function (x) {
  inherits(x, 'transition')
}


stage_matrixCheck <- function (x) {
  stopifnot(ncol(x) == nrow(x))
  stopifnot(is.matrix(x))
  stopifnot(all(is.finite(x)))
}

transitionCheck <- function (x) {
  stopifnot(inherits(x,'transition'))
  stopifnot(is.list(x))
}

subTransition <- function (x, i) {
  attrib <- attributes(x)
  attrib$names <- attrib$names[i]
  x <- x[i]
  attributes(x) <- attrib
  return (x)
}

dots <- function(...) {
  eval(substitute(alist(...)))
}
