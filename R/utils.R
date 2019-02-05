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

set.seed(123)
population <- runif(10) + runif(10)

round_pop <- function (population) {

  population_min <- floor(population)

  if (steps_stash$demo_stochasticity == "full") {
    return(stats::rmultinom(1, sum(population), population)[, 1])
  }
  
  if (steps_stash$demo_stochasticity == "deterministic_redistribution") {
    n <- length(population)
    k <- round(sum(population)) - sum(population_min)
    cutoff <- sort(population, partial = n - k)[n - k]
    idx <- which(population > cutoff)
    population_min[idx] <- population_min[idx] + 1
    return(population_min)
  }
  
  if (steps_stash$demo_stochasticity == "stochastic_redistribution") {
    population_extra <- population - population_min
    population_extra[] <- stats::rbinom(length(population_extra), 1, population_extra[])
    return(population_min + population_extra)
  }
  
  if (steps_stash$demo_stochasticity == "none") return(population)
  
}