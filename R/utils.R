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

round_pop <- function (population) {
  population_min <- floor(population)
  population_extra <- population - population_min
  population_extra[] <- stats::rbinom(length(population_extra), 1, population_extra[])
  population <- population_min + population_extra
  population
}