#' @title dispersal class for metapoplation
#' @rdname dispersal 
#' @name dispersal
#' @description creates a function that governs dispersal capacity of stages in a population. The input is a list which contains 
#'  the first list is a dispersal kernal \code{alpha}  value for each stage, the second named list \code{probability} 
#' is the proportion of that stage will disperse. For example a probability of 0.2 for stage larvae
#' means a random 20% of larve will a try and disperse to patches, the distance they can disperse is governed 
#' by the dispersal kernal (alpha). If params = NULL, a dispersal kernal of 1 is given to all stages, 
#' and all stages will attempt to undertake dispersal.
#' @param params named list w
## and disp_fun a characture 'H' uses hanski(1994), if 'S' uses shaw(1995).
#' @export
#' @examples 
#' params <- list(alpha=list('larval'=2,'juvenile'=0,'adult'=3),
#'                probability=list('larval'=0.2,'juvenile'=0,'adult'=0.6))  
#'                
#' dispersal(params)
#' 
#' dispersal(NULL)

dispersal <- function(params,...){
  switch(class(params[1]),
         NULL = dispersalDefault(...),
         list = dispersal_core(params,...))
}

#' @rdname dispersal
#' @name d
#' @export
#' @examples 
#' d(NULL)
#' 
d <- dispersal

dispersal_core <- function(params,...)
  {
  stopifnot(is.list(params))
  ds <- function(habitat) {
          disp <- lapply(params$alpha,function(x,dist)exp(-x*dist),dist=habitat$distance)
          disp <- lapply(disp,function(x)sweep(x,1,rowSums(x),'/'))
          disp
  }
  params$disp <- ds
  ds <- as.dispersal(params)
  return(ds)
}


dispersalDefault <- function (...) {
  params_list <- list(alpha = list("alpha"=1))
  d <- dispersal_core(params_list)
  return (d)
}

#' @rdname dispersal
#' @param x an object to be tested as a dispersal transfun object
#' @export
#' @examples
#' is.dispersal(disp)
is.dispersal <- function (x) inherits(x, 'dispersal')

#' @rdname dispersal
#' @param x an object to be tested as a dispersal transfun object
#' @export
as.dispersal <- function (x) {
  if (!is.transition(x)) {
      class(x) <- c('transition', class(x))
    }
    if (!is.dispersal(x)) {
      class(x) <- c('dispersal', class(x))
    }
    return (x)
}

#' @rdname dispersal
#' @param x an object to print or test as a transition object
#' @param \dots further arguments passed to or from other methods.
#' @export
#' @examples
#' # print method
#' dispersal(params)
#' print(d(NULL))
#'
print.dispersal <- function(x,...){
  disp_info <- list()
  for (i in base::seq_along(x)[-length(x)]) {
  disp_info[[i]]  <- base::paste(
      base::names(x[[i]]),
      base::format(x[[i]], scientific = FALSE),
      sep = " = ",
      collapse = "\n "
    ) 
  };
  if(length(disp_info)==2) text <- sprintf('dispersal function with disperal kernal of:\n %s\n and probability of dispersal of:\n %s\n for stages.',disp_info[[1]],disp_info[[2]])
  else text <- sprintf('dispersal function with disperal kernal of:\n %s\n for stages.',disp_info[[1]])
  cat(text)
}

probdisp <- function (disp, habitat) {
  # get expected dispersal fraction from a probability and a dispersal transfun.
  # dispersal should have diagonal giving probability of staying, off-diagonals
  # giving probability of moving to each other patch, w/ all rows summing to 1.
  # probability is probability of leaving (1-probability of staying).

  # work out which way round
  disps <- disp$disp(habitat)
  probs <- disp$probability

  # multiply each row by the dispersal probability
  disps <- mapply(FUN=function(x,y) sweep(x,1,y,'*'),x=disps, y=probs,SIMPLIFY = FALSE)

  # add fraction not attempting dispersal back onto diagonal
  for(i in seq_along(disps))diag(disps[[i]]) <- diag(disps[[i]]) + 1 - probs[[i]]
  return (disps)
}