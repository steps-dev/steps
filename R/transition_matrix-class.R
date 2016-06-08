#' @title transition matrix objects 
#' @name as.transition_matrix
#' @rdname transition_matrix
#' @param x For as.transition_matrix, x is a square matrix, that has transition states between population stages.
#' @param names.st string of names for each stage
#' @param ... other function calls.
#' @return An object of class tmatrix, i.e, resulting from as.transition_matrix.
#' @author Skipton Woolley
#' @export
#' @examples 
#' mat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' tmat <- as.transition_matrix(mat)

as.transition_matrix <- function(x, names.st=NULL,...){
    x <- base::as.matrix(x,...)
    if(base::diff(base::dim(x)) !=0) stop("Needs to be a square matrix with transition probabilities between each stage.")
    di <- base::dim(x)[[1]]
    m.names <- base::dimnames(x)[[1]]
    if(base::is.null(m.names)) m.names <- names.st
    if(base::is.null(m.names)) m.names <- base::paste0("stage.",1:di)
    base::dimnames(x) <- base::list(m.names, m.names)
    base::class(x)<-c("transition_matrix", class(x))
    return(x)
}

#' @rdname transition_matrix
#' @name is.transition_matrix
#' @export
#' @examples
#' is.transition_matrix(tmat)
is.transition_matrix <- function (x) {
  inherits(x, 'transition_matrix')
}

#' @rdname transition_matrix
#' @param object an object of \code{transition_matrix} class
#' @export
#' @description prints the main parameters of the transition matrix: the finite rate of increase ("lambda"), the stable stage distribution, the reproductive value and the sensitivities and elasticities matrices.
#' @examples
#' mat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' tmat <- as.transition_matrix(mat)
#' summary(tmat)
summary.transition_matrix <-
  function(object,...){
    name.mat<-deparse(substitute(object))
    x <- object
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
    class(result)=c("summary.transition_matrix", class(result))
    return (result)
  }


#' @rdname transition_matrix
#' @importFrom igraph graph.adjacency
#' @export
#' @author Nick Golding
#' @examples 
#' mat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' tmat <- as.transition_matrix(mat)
#' plot(tmat)

plot.transition_matrix <- function (x, ...) {
  # plot a dynamic using igraph
  
  # extract the transition matrix & create an igraph graph x
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
