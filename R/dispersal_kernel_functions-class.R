#' Create a dispersal function
#'
#' @description A dispersal kernal function is mathematical representation of how species redistribute
#' across the landscape.
#'  
#' More common kernels are built-in for the user to select (see exponential_dispersal_kernel), however,
#' a user may also provide a custom written dispersal kernel.
#' 
#' @rdname dispersal_function
#'
#' @param distance_decay a parameter to control the rate at which the population disperses with distance
#' @param normalize should the normalising constant be used - default is FALSE.
#' 
#' @return An object of class \code{dispersal_function}
#' 
#' @export
#'
#' @examples
#' 
#' test_dispersal_function <- exponential_dispersal_kernel()

exponential_dispersal_kernel <- function (distance_decay = 0.5, normalize = FALSE) {
  if (normalize) {
    fun <- function (r) (1 / (2 * pi * distance_decay ^ 2)) * exp(-r/distance_decay)
  } else {
    fun <- function (r) exp(-r/distance_decay)
  }
  as.dispersal_function(fun, info = list(distance_decay = distance_decay,
                                         normalize = normalize))
}

# #' @rdname dispersal_function
# #'
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

as.dispersal_function <- function (dispersal_function, info = NULL) {
  as_class(dispersal_function, "dispersal_function", "function", info = info)
}

print.dispersal_function <- function (x, ...) {
  print_info(x)
}
