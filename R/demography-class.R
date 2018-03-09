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

build_demography <- function (transition_matrix, dispersal_parameters) {
  demography <- list(transition_matrix = transition_matrix,
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
#' @examples

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

##########################
### internal functions ###
##########################

demographyCheck <- function (x) {
  stopifnot(inherits(x,'demography'))
  stopifnot(is.list(x))
}
