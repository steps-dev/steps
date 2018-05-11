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

# given a 3D array, and a vector of values (or cell numbers if index = TRUE) in
# the first slice (a matrix), return a vector of the corresponding
# vcalues/indices for corresponding elements in the whole array
replicate_values <- function (values, array, index = FALSE) {
  
  # get relevant dimensions fomr the 3D array
  dims <- dim(array)
  ndim <- length(dims)
  
  if (ndim == 3) {
    
    nelem <- prod(dims[1:2])
    replicated <- dims[3]
    
    # handle differently if we're replicating an index
    if (index) {
      
      # how much to pad the index to get corresponding elements in each slice  
      addition <- nelem * (seq_len(replicated)  -1)
      # get indices for the whole array
      values <- as.numeric(outer(values, addition, FUN = "+"))
      
    } else {
      values <- rep(values, replicated)
    }
    
  } else if (ndim != 2) {
    stop ("Check dimensions of input object")
  }
  
  values
  
}
