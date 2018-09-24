#' Create a demography object to use in a state object.
#'
#' A demography object contains information on how population
#' demographics change in space and time.
#' 
#' A demography object is a sub-component of a \link[steps]{state} object
#' and is modified in each timestep of a simulation. During a simulation,
#' a demography object tracks changes in life-stage matrices (either global
#' for the entire landscape or local to grid cells) based on demographic
#' dynamic functions selected or created by the user. 
#' 
#' @rdname demography
#' 
#' @param transition_matrix A symmetrical age-based (Leslie) or stage-based
#' population structure matrix.
#' @param scale The scale to which the transition matrix is applied - either
#' 'global' for a landscape-wide application or 'local' for a grid cell-based
#' applications. The default value is 'global'.
#' @param habitat_suitability A spatial raster (grid cell-based) layer that
#' is of the appropriate extent and resolution for the simulation (required
#' if 'local' scale is specified).
#' @param misc Miscellaneous inputs used to modify the demography object in a
#'  simulation. Note, this is where the user may store objects that are called
#'  upon and modified by demographic dynamic functions (also created by the user).
#' @param x A demography object (test or print).
#' @param object A demography object (summary).
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{demography}
#' 
#' @export
#'
#' @examples
#' 
#' library(steps)
#' library(raster)
#' 
#' # Construct the demography object
#' dem <- demography(transition_matrix = egk_mat)

demography <- function (transition_matrix, 
                        scale = 'global', 
                        habitat_suitability = NULL, 
                        misc = NULL, 
                        ...) {
  x <- transition_matrix
  stage_matrixCheck(x)
  di <- base::dim(x)[[1]]
  m.names <- base::dimnames(x)[[1]]
  if(base::is.null(m.names)) m.names <- base::paste0("stage.",1:di)
  base::dimnames(x) <- base::list(m.names, m.names)
  
  if ( scale == 'local' ) {
    if ( !inherits(habitat_suitability,c("RasterLayer")) ) {
      stop("A raster layer must be specified if storing local (cell-based) transition matrices")
    }
    
    ncells <- length(habitat_suitability)
    
    # set up a matrix of (ncell)(rows)*(nstage*nstage)(cols) and fill with original stage matrix
    #all_stage_matrices <- matrix(rep(c(x),ncells),ncells,length(c(x)),byrow = TRUE)
    #colnames(all_stage_matrices)<-paste0(rep(m.names,each=di),1:di)
    
    # setup an array of matrices that is nstage*nstage*ncells:
    all_stage_matrices <- array(NA, dim = c(di, di, ncells))
    all_stage_matrices[ , , seq_len(ncells)] <- transition_matrix
    
    demography <- list(global_transition_matrix = x,
                       local_transition_matrix = all_stage_matrices,
                       misc = misc)
  }else{
    demography <- list(global_transition_matrix = x,
                       local_transition_matrix = NULL,
                       misc = misc)    
  }
  
  as.demography(demography)
}

#' @rdname demography
#'
#' @export
#' 
#' @examples
#' 
#' # Test if object is of the type 'demography'
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
#' print(test_demography)

print.demography <- function (x, ...) {
  if (!is.null(x$local_transition_matrix)) {
    
    cat("This is a demography object that contains", length(x$local_transition_matrix),
        " independent transition matrices - one for each grid cell in the landscape.")
    
  } else {
    
    cat("This is a demography object that contains a single transition matrix",
        " that will be applied similarly to all grid cells across the landscape.")
    
  }
}

#' @rdname demography
#'
#' @export
#'
#' @examples
#' 
#' # Print a summary of 'demography' object attributes
#' summary(test_demography)

summary.demography <- function (object, ...) {

  x <- object$global_transition_matrix
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

as.demography <- function (demography) {
  as_class(demography, "demography", "list")
}

stage_matrixCheck <- function (x) {
  if (!is.matrix(x)) {
    stop("A matrix object is required.")
  }
  
  if (ncol(x) != nrow(x)) {
    stop("A square matrix with stage probabilities between each stage is required.")
  }
  
  if (!all(is.finite(x))) {
    stop("All values in matrix are required to be either zero or positive and finite.")
  }
}
