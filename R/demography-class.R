#' Create a demography object to define population transition matrices
#'
#' @param transition_matrix A symmetrical stage-based population structure matrix
#' @param dispersal_parameters Specifications for dispersal in the landscape
#'
#' @return An object of class \code{demography}
#' @export
#'
#' @examples
#' 
#' library(raster)
#' library(dhmpr)
#' 
#' test_demography <- build_demography(fake_transition_matrix(4), rlnorm(1))

build_demography <- function (transition_matrix, transition_matrix_sd=0, dispersal_parameters) {
  x <- transition_matrix
  if(base::diff(base::dim(x)) !=0) stop("Needs to be a square matrix with stage probabilities between each stage.")
  stage_matrixCheck(x)
  di <- base::dim(x)[[1]]
  m.names <- base::dimnames(x)[[1]]
  if(base::is.null(m.names)) m.names <- base::paste0("stage.",1:di)
  base::dimnames(x) <- base::list(m.names, m.names)
  demography <- list(transition_matrix = transition_matrix,
                     transition_matrix_sd = transition_matrix_sd,
                     dispersal_parameters = dispersal_parameters)
  set_class(demography, "demography")
}

#' Print details of a demography object
#'
#' @param x an object to print or test as an demography object
#' @param ... further arguments passed to or from other methods
#'
#' @export
#'
# @examples
# test_demography <- build_demography(fake_transition_matrix(4), rlnorm(1))
# print(test_demography)

print.demography <- function (x, ...) {
  cat("This is a demography object")
}

#' Show details of demography object
#'
#' @param x a demography object to summarise
#' @param ... further arguments passed to or from other methods 
#'
#' @export
#'
##' @examples
##' summary(demo)

summary.demography <- function (x,...) {
  demographyCheck(x)
  x <- x$transition_matrix
  di <- base::dim(x)[1]
  ea <- base::eigen(x)
  lambda <- base::abs(ea$values[1]) 
  ssd <- base::abs(ea$vectors[,1]/base::sum(ea$vectors[,1]) ) 
  ae <- base::eigen(base::t(x))
  vr <- base::abs(ae$vectors[,1]/ae$vectors[1,1] )
  sensitivity <- (vr  %*%  base::t(ssd))  / (base::t(vr) %*% ssd)[1,1]
  elasticity <- sensitivity * x / lambda
  
  #### Add stochasticity summary...
  result<- list(lambda = lambda,
                stable.stage.distribution = ssd,
                reproductive.value = vr,
                sensitivity = sensitivity,
                elasticity = elasticity
                )
  
  return (result)
}

#' Plot details of demography object
#'
#' @param x a demography object to summarise
#' @param ... further arguments passed to or from other methods 

#' @importFrom igraph graph.adjacency
#' @export
#' 
#' @examples 
##' plot(demo)

plot.demography <- function (x, ...) {
  # plot a dynamic using igraph
  
  # extract the stage matrix & create an igraph graph x
  graphics::par(mar=c(2,4,4,2))
  demographyCheck(x)
  x <- x$global_stage_matrix
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

#' Verify demography object
#'
#' @param x a demography object to summarise
#'
#' @export
#' 
#' @examples
##' is.demography(demo)

is.demography <- function (x) {
  inherits(x, 'demography')
}

##########################
### internal functions ###
##########################

stage_matrixCheck <- function (x) {
  stopifnot(ncol(x) == nrow(x))
  stopifnot(is.matrix(x))
  stopifnot(all(is.finite(x)))
}

demographyCheck <- function (x) {
  stopifnot(inherits(x,'demography'))
  stopifnot(is.list(x))
}
