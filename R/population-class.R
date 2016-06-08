#' population class
#' @title population objects
#' @name population
#' @rdname population
#' @param x a vector, data.frame or matrix of population(s) numbers for each stage(s) and each patch(s)
#' @export
#' @examples 
#' # starting population numbers for each step in the demographic function
#' population <- as.population(c(80,30,10))
#' population <- as.population(t(rmultinom(10, size = 100, prob = c(0.8,0.2,0.01))))

as.population <- function(x,...){
  if(base::is.null(base::dim(x))){
    names(x) <- base::paste0("stage",seq_along(x))
  } else { 
    x <- base::as.matrix(x,...)
  ns <- base::dim(x)[[2]]
  np <- base::dim(x)[[1]]
  s.names <- base::dimnames(x)[[2]]
  p.names <- base::dimnames(x)[[1]]
  if(base::is.null(s.names)) s.names <- base::paste0("stage",1:ns)
  if(base::is.null(p.names)) p.names <- base::paste0("patch",1:np)
  base::dimnames(x) <- base::list(p.names, s.names)
  }
  base::class(x)<-c("population", class(x))
  return(x)
}

#' @rdname population
#' @name is.population
#' @export
#' @examples
#' population <- as.population(c(80,30,10))
#' is.population(population)
is.population <- function (x) {
  inherits(x, 'population')
}

