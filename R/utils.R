as_class <- function (object, name, type = c("function", "list")) {
  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))
  invisible(object)
}

round_pop <- function (population) {
  
  population_min <- floor(population)
  
  if (steps_stash$demo_stochasticity == "full") {
    if (sum(population) == 0) return(population)
    return(rmultinom_large_int(population)[, 1])
  }

  if (steps_stash$demo_stochasticity == "none") return(population)
  
}

get_carrying_capacity <- function (landscape, timestep) {
  
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
    
    names(res) <- paste0("Carrying_Capacity_", timestep)
    
  } else {
    
    # otherwise, we don't support it
    stop ("invalid carrying capacity argument",
          call. = FALSE)
    
  }
  
  res
  
}

not_missing <- function (raster) {
  which(!is.na(raster::getValues(raster)))
}

warn_once <- function (condition, message, warning_name) {
  if (condition & !isTRUE(steps_stash[[warning_name]])) {
    warning(message, call. = FALSE)
    steps_stash[[warning_name]] <- TRUE
  }
}

rmultinom_large_int <- function (population) {
  
  total <- round(sum(population))
  
  if (total > .Machine$integer.max) {
    times <- total %/% .Machine$integer.max
    extra <- total %% .Machine$integer.max
    pop <- rep(0, length(population))
    
    for (i in seq_len(times)) {
      pop <- pop + stats::rmultinom(1, .Machine$integer.max, population)
    }
    
    pop <- pop + stats::rmultinom(1, extra, population)
    
  } else {
    
    pop <- stats::rmultinom(1, total, population)
    
  }
  
  pop
  
}

pretty_int <- function (...) {
  at <- pretty(...)
  at <- at[at %% 1 == 0]
  at[at != 0]
}

int_or_proper_length_vector <- function (input, n_stages, parameter) {
  if (length(input) != 1 & length(input) != n_stages) {
    stop(paste0("Please provide either a single number or vector of ",
                "numbers that matches the number of life-stages in the ",
                parameter,
                " parameter."))
  }
  if (length(input) == 1) {
    input <- rep(input, n_stages)
  }
  input
}

global_object_error <- function(error) {
  
  # see if there's a missing object
  something_missing <- grepl("could not find", error$message)
  
  if (something_missing) {
    message <- paste(error$message,
                     "\n\nit looks like the future package can't find an object or",
                     "function you are using, you can pass it in via the",
                     "future.globals argument to steps::simulation")
  } else {
    message <- error$message
  }
  
  stop (message, call. = FALSE)
}