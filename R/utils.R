#' @importFrom magrittr %>% %<>%

# to use magrittr shortcut
utils::globalVariables(".")

print_info <- function (x) {
  info <- attr(x, "info")
  if (!is.null(info) && !identical(info, list())) {
    for (member in names(info)) {
      obj <- info[[member]]
      if (!is.null(obj)) {
        if("function" %in% class(obj)) msg <- paste0(toupper(member), " FUNCTION:\n\n")
        else msg <- paste0(member, ":\n\n")
        cat(msg)
        print(obj)
        cat("\n")
      }
    }
  }
}

as_class <- function (object, name, type = c("function", "list"), info = list()) {
  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))
  attr(object, "info") <- info
  object
}

round_pop <- function (population) {

  population_min <- floor(population)

  if (steps_stash$demo_stochasticity == "full") {
    if (sum(population) == 0) return(population)
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

get_carrying_capacity <- function(landscape, timestep) {
  
  cc <- landscape$carrying_capacity
  if (is.null(cc)) {
    
    # if there's no carrying capacity specified, return a NULL
    res <- NULL
    
  } else if (inherits(cc, "RasterLayer")) {
    
    
    # if there's a carrying capacity raster, use that
    
    # # in a future set up where lots of carrying capacity rasters could be passed in
    # if (raster::nlayers(cc) > 1) {
    #   res <- cc[[timestep]])
    # } else {
    #   res <- cc
    # }
    
    res <- cc
    
  } else if (is.function(cc)) {
    
    # if it's a function, run it on landscape
    res <- cc(landscape, timestep)
    
  } else {
    
    # otherwise, we don't support it
    stop ("invalid carrying capacity argument",
          call. = FALSE)
    
  }
  
  res
  
}

not_missing <- function(raster) {
  which(!is.na(raster::getValues(raster)))
}
