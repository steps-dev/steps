#' @title event class for dynamic metapopulation class models.
#' @rdname event 
#' @description This class is a way to order, schedule and run events which will have
#' an effect on the dynamic metapopulation models. Events can be temporally assigned. 
#' But need to be ordered with respect to metapopulation steps.
#' Events can have different temporal requirements. For example, management and demographic events
#' might be on a yearly scale, were fire might be weekly.

#' @export
#' @rdname event 
order_events <- function(){
  
}

#' @export
#' @rdname event 
schedule_event <- function(){
  
}

#' @export
#' @rdname event 
do_event <- function(){
  
}
#' Event order
#' Assign event order: 1 = first (highest); 5 = normal; 10 = last (lowest).
#'
#' @return A numeric.
#' @export
#' @rdname event
#'
.first <- function() {
  .highest()
}

#' @export
#' @rdname event
.highest <- function() {
  return(1)
}

#' @export
#' @rdname event
.last <- function() {
  .lowest()
}

#' @export
#' @rdname event
.lowest <- function() {
  return(10)
}

#' @export
#' @rdname event
.normal <- function() {
  5
}

