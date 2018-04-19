#' @importFrom magrittr %>% %<>%

# to use magrittr shortcut
utils::globalVariables(".")

# to set object classes
# set_class <- function (x, class) {
#   class(x) <- c(class, class(x))
#   x
# }

as_class <- function (object, name, type = c("function", "list")) {
  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))
  object
}

# create a transition matrix for testing
fake_transition_matrix <- function (n_stages) {
  
  survival <- stats::runif(n_stages, 0.7, 0.9)
  growth <- stats::runif(n_stages - 1, 0.5, 0.7)
  recruitment <- stats::rlnorm(1)
  
  # base matrix 
  transition_matrix <- diag(n_stages)
  
  # add growth
  growth_idx <- which(row(transition_matrix) == col(transition_matrix) + 1)
  transition_matrix[growth_idx] <- growth
  
  # columns sum to 1
  sums <- colSums(transition_matrix)
  transition_matrix <- sweep(transition_matrix, 2, sums, FUN = "/")
  
  # apply survival
  transition_matrix <- sweep(transition_matrix, 2, survival, FUN = "*")
  
  # add recruitment
  transition_matrix[1, n_stages] <- recruitment
  
  transition_matrix
}
