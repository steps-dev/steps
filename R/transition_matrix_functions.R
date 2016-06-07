#' Transition population matrix functions

#' @title transition matrix objects 
#' @name trans_matrix
#' @rdname trans_matrix
#' @param x	For as.trans_matrix, x is a square matrix, that has transition states between population stages.
#' @return An object of class tmatrix, i.e, resulting from as.trans_matrix.
#' @author Skipton Woolley
#' @export
#' @examples 
#' mat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' tmat <- as.trans_matrix(mat)

as.trans_matrix <- function(x, names.st=NULL,...){
    x <- base::as.matrix(x,...)
    if(base::diff(base::dim(x)) !=0) stop("Needs to be a square matrix with transition probabilities between each stage.")
    di <- base::dim(x)[[1]]
    m.names <- base::dimnames(x)[[1]]
    if(base::is.null(m.names)) m.names <- names.st
    if(base::is.null(m.names)) m.names <- base::paste0("stage.",1:di)
    base::dimnames(x) <- base::list(m.names, m.names)
    base::class(x)<-c("trans_matrix", class(x))
    return(x)
}

#' @rdname trans_matrix
#' @param x a trans_matrix object
#' @export
#' @examples
#' is.trans_matrix(tmat)
is.trans_matrix <- function (x) {
  inherits(x, 'trans_matrix')
}

#' #' @rdname trans_matrix
#' #' @param x a trans_matrix object
#' #' @param ... other function calls from igraph.
#' #' @importFrom igraph graph.adjacency
#' #' @export
#' get_trans_matrix_graph <-function(x,...)
#' {
#'   tmg <- igraph::graph.adjacency(adjmatrix=x, weighted=TRUE, mode="directed",...)
#'   return(tmg)
#' }

