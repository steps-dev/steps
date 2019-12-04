#' Create a dispersal function
#'
#' A dispersal kernel function is a mathematical representation of how species redistribute
#' across the landscape.
#'  
#' A common dispersal kernel is provided in the software for the user to select, however,
#' a user may also provide a custom written dispersal kernel.
#' 
#' @name dispersal_kernel
#' @seealso
#' \itemize{
#'   \item{\link[steps]{exponential_dispersal_kernel}) for a (negative) exponential dispersal
#'   kernel}
#'   }
NULL

#' Negative exponential dispersal kernel
#'
#' This function determines the proportion of redistribution based on distance.
#' 
#' @param distance_decay (exponential dispersal parameter) controls the rate at which the population disperses with distance
#' @param normalize (exponential dispersal parameter) should the normalising constant be used - default is FALSE.
#' 
#' @return An object of class \code{dispersal_function}
#' 
#' @export
#'
#' @examples
#'
#' \dontrun{
#' dists <- seq(0, 100, 1)
#' 
#' exp_dispersal_fun <- exponential_dispersal_kernel(distance_decay = 50)
#' 
#' plot(dists, exp_dispersal_fun(dists), type = 'l')
#' }

exponential_dispersal_kernel <- function (distance_decay = 0.5, normalize = FALSE) {
  if (normalize) {
    fun <- function (r) (1 / (2 * pi * distance_decay ^ 2)) * exp(-r / distance_decay)
  } else {
    fun <- function (r) exp(-r / distance_decay)
  }
  as.dispersal_function(fun)
}

# #' @param x an object to print or test as a dispersal_function object
# #' @param ... further arguments passed to or from other methods
# #'
# #' @export
# #'
# #' @examples
# #'
# #' print(test_dispersal_function)
# 
# print.dispersal_function <- function (x, ...) {
#   cat("This is a dispersal_function object")
# }

##########################
### internal functions ###
##########################

as.dispersal_function <- function (dispersal_function) {
  as_class(dispersal_function, "dispersal_function", "function")
}
