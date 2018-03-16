#' Create a demography object to use in a state object
#'
#' @description A demography object is used to store information on how populations change in space and time.
#' This includes life-stage matrices and parameters to control dispersal.
#' It is a sub-component of a \link[dhmpr]{state} object and is modified in each timestep of an experiment.
#' 
#' @rdname demography
#' 
#' @param transition_matrix a symmetrical stage-based population structure matrix
#' @param type scale of transition matrix - either 'global' or 'local' (cell-based)
#' @param habitat_suitability habitat suitability raster layer (required if 'local' type is specified)
#' @param dispersal_parameters specifications for dispersal in the landscape
#' @param misc miscellaneous inputs used to modify demography
#' @param x a demography object
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of class \code{demography}
#' 
#' @export
#'
#' @examples
#' 
#' library(dhmpr)
#' library(raster)
#' 
#' # Use a built-in function to generate a four life-stage transition matrix
#' 
#' mat <- fake_transition_matrix(4)
#' 
#' # Provide a list of dispersal parameters
#' 
#' params <- list(rlnorm(1),c(1:4))
#' 
#' # Construct the demography object.
#' 
#' test_demography <- build_demography(transition_matrix = mat, dispersal_parameters = params)

build_demography <- function (transition_matrix, type = 'global', habitat_suitability = NULL, dispersal_parameters, misc = NULL, ...) {
  x <- transition_matrix
  stage_matrixCheck(x)
  di <- base::dim(x)[[1]]
  m.names <- base::dimnames(x)[[1]]
  if(base::is.null(m.names)) m.names <- base::paste0("stage.",1:di)
  base::dimnames(x) <- base::list(m.names, m.names)
  
  if (type == 'local' ) {
    if (!inherits(habitat_suitability,c("RasterLayer"))) stop("A raster layer must be specified if storing local (cell-based) transition matrices")
    ncells <- length(habitat_suitability)
    
    # set up a matrix of (ncell)(rows)*(nstage*nstage)(cols) and fill with original stage matrix
    all_stage_matrices <- matrix(rep(c(x),ncells),ncells,length(c(x)),byrow = TRUE)
    colnames(all_stage_matrices)<-paste0(rep(m.names,each=di),1:di)
    
    demography <- list(global_transition_matrix = x,
                       local_transition_matrix = all_stage_matrices,
                       dispersal_parameters = dispersal_parameters,
                       misc = misc)
  }else{
    demography <- list(global_transition_matrix = x,
                       dispersal_parameters = dispersal_parameters,
                       misc = misc)    
  }

  set_class(demography, "demography")
}

#' @rdname demography
#'
#' @export
#' 
#' @examples
#' 
#' # Test if object is of the type 'population'
#' 
#' is.demography(test_demography)

is.demography <- function (x) {
  inherits(x, 'demography')
}

#' @rdname demography
#'
#' @export
#'
#' @examples
#' 
#' # Print information about the 'demography' object
#' 
#' print(test_demography)

print.demography <- function (x, ...) {
  cat("This is a demography object")
}

#' @rdname demography
#'
#' @export
#'
#' @examples
#' 
#' # Print a summary of 'demography' object attributes
#' 
#' summary(test_demography)

summary.demography <- function (x,...) {

  x <- x$global_transition_matrix
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

#' @rdname demography
#'
#' @importFrom igraph graph.adjacency
#' 
#' @export
#' 
#' @examples
#' 
#' # Plot the 'demography' object
#' 
#' plot(test_demography)

plot.demography <- function (x, ...) {
  # plot a dynamic using igraph
  
  # extract the stage matrix & create an igraph graph x
  graphics::par(mar=c(2,4,4,2))
  x <- x$global_transition_matrix
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

##########################
### internal functions ###
##########################

stage_matrixCheck <- function (x) {
  if (!is.matrix(x)) stop("A matrix object is required")
  if (ncol(x) != nrow(x)) stop("A square matrix with stage probabilities between each stage is required")
  if (!all(is.finite(x))) stop("All values in matrix are required to be either zero or positive and finite")
}
