#' @title fire spread object for dynamic habitat meta-populations.
#' @name fire_spread
#' @rdname fire_spread
#' @description This is an example function which we can import as a module to manipulate the habitat. 
#' This function uses cellular automata to spread fire across the habitat. These burnt cells can be used to reduce (or increase)
#' the probability of occurrence within the \code{habitat} object and thus re-shape or change patch suitability. 
#' Which will have an effect on the metapopulation models.
#' @export
#'
fire_spread <- function(habitat,
                        fire_start_location = NA_real_,
                        prob = 0.23, # could replace this with cell probabilities.
                        continue_to_burn_prob = 0,
                        max_cells = 1e8L,
                        directions = 8L,
                        iterations = 1e6L,
                        return_table = FALSE,
                        id = FALSE,
                        spread_prob_after_ignition = NA_real_,
                        fire_history = NA) {
    
    habitat <- suitability(habitat)
    fire_historyExists <- is(fire_history, "data.table")
    if (!is(spread_prob_after_ignition, "Raster")) {
      if (is.na(spread_prob_after_ignition)) {
        spread_prob_after_ignition <- prob
      }
    }
    ### should sanity check map extents
    if (any(is.na(fire_start_location)))  {
      # start it in the centre cell, if there is no fire_history
      if (!fire_historyExists)
        fire_start_location <- (nrow(habitat)/2L + 0.5) * ncol(habitat)
    }
    
    if (fire_historyExists) {
      fire_start_location <- fire_start_location[!(fire_start_location %in% fire_history[,indices])] # keep these for later
      initialfire_start_location <- fire_start_location
    } else {
      initialfire_start_location <- fire_start_location
    }
    
    spreads <- vector("integer", ncell(habitat))

    n <- 1L
    
    if (id | return_table) { # give values to spreads vector at initialfire_start_location
        spreads[fire_start_location] <- 1L:length(fire_start_location)
      } else {
        spreads[fire_start_location] <- n
      }
    
    # Convert NAs to 0 on the spread_prob_after_ignition Raster
    if (is(spread_prob_after_ignition, "Raster")) {
      # convert NA to 0s
      spread_prob_after_ignition[is.na(spread_prob_after_ignition)] <- 0L
    } 
    
    # Convert NAs to 0 on the prob Raster
    if (is(prob, "Raster")) {
      # convert NA to 0s
      prob[is.na(prob)] <- 0L
    } 
    
    if (fire_historyExists) {
      if (sum(colnames(fire_history) %in% c("indices", "id", "active", "initialLocus")) == 4) {
          spreads[fire_start_location] <- spreads[fire_start_location] + fire_history[, max(id)] # reassign old ones
          spreads[fire_history[,indices]] <- fire_history[, id]
          fire_start_location <- c(fire_history[active == TRUE, indices], fire_start_location) %>% na.omit()
        } else {
          stop("fire_history must have at least columns: ",
               "indices, id, active, and initialLocus.")
        }
      }
    
    if (any(fire_start_location > ncell(habitat))) stop("fire_start_location indices are not on habitat")
    
    ## Recycling max_cells as needed
    if (any(!is.na(max_cells))) {
      if (!is.integer(max_cells)) max_cells <- floor(max_cells)
      if (fire_historyExists) {
        sizeAll <- fire_history[, list(len = length(initialLocus)), by = id]
        max_cells <- rep_len(max_cells, length(initialfire_start_location) + NROW(sizeAll))
        size <- c(sizeAll[, len], rep_len(1L, length(initialfire_start_location)))
      } else {
        max_cells <- rep_len(max_cells, length(fire_start_location))
        size <- rep_len(1L, length(fire_start_location))
      }
    } else {
      max_cells <- ncell(habitat)
      size <- length(fire_start_location)
    }
    
    # while there are active cells
    while (length(fire_start_location) & (n <= iterations) ) {
      
      # identify neighbours
      if (id | return_table) {
          potentials <- adj(habitat, fire_start_location, directions, pairs = TRUE)
        } else {
          # must pad the first column of potentials
          potentials <- cbind(NA, adj(habitat, fire_start_location, directions, pairs = FALSE))
      }
      
      
      # keep only neighbours that have not been spread to yet
      potentials <- potentials[spreads[potentials[, 2L]] == 0L, , drop = FALSE]

      if (n == 2) {
        prob <- spread_prob_after_ignition
      }
      
      if (is.numeric(prob)) {
        if (n == 1 & fire_historyExists) { # need cell specific values
          probs <- rep(prob, NROW(potentials))
          prevIndices <- potentials[, 1L] %in% fire_history[active == TRUE, indices]
          probs[prevIndices] <- spread_prob_after_ignition
        } else {
          probs <- prob
        }
      } else {
        if (n == 1 & fire_historyExists) { # need cell specific values
          probs <- prob[potentials[, 2L]]
          prevIndices <- potentials[, 1L] %in% fire_history[active == TRUE, indices]
          probs[prevIndices] <- spread_prob_after_ignition
        } else {
          probs <- prob[potentials[, 2L]]
        }
      }

      if (any(probs < 1)) {
        potentials <- potentials[runif(NROW(potentials)) <= probs, , drop = FALSE]
      }
      potentials <- potentials[sample.int(NROW(potentials)), , drop = FALSE] # random ordering so not always same
      potentials <- potentials[!duplicated(potentials[, 2L]), , drop = FALSE]

      
      if (length(potentials) > 0) {# potentials can become zero because all active cells are edge cells
        
        events <- potentials[, 2L]
        
        # Implement max_cells
        if (length(max_cells) == 1L) {
          len <- length(events)
          if ((size + len) > max_cells) {
            keep <- len - ((size + len) - max_cells)
            samples <- sample(len, keep)
            events <- events[samples]
            potentials <- potentials[samples, , drop = FALSE]
          }
          size <- size + length(events)
        } else {
          len <- tabulate(spreads[potentials[, 1L]], length(max_cells))
          if ( any( (size + len) > max_cells & size <= max_cells) ) {
            whichID <- which(size + len > max_cells)
            toRm <- (size + len)[whichID] - max_cells[whichID]
            for (i in 1:length(whichID)) {
              thisID <- which(spreads[potentials[, 1L]] == whichID[i])
              potentials <- potentials[-sample(thisID, toRm[i]), , drop = FALSE]
            }
            events <- potentials[, 2L]
          }
          size <- pmin(size + len, max_cells) ## Quick? and dirty. fast but loose (too flexible)
        }
        
        # increment iteration
        n <- n + 1L
        
        if (length(events) > 0) { # place new value at new cells that became active
            if (id | return_table) {
              spreads[events] <- spreads[potentials[, 1L]]
            } else {
              spreads[events] <- n
            }
          }
        
        # remove the cells from "events" that push it over max_cells
        if (length(max_cells) > 1L) {
          if (exists("whichID", inherits = FALSE)) {
              max_cellsKeep <- !spreads[events] %in% whichID
              events <- events[max_cellsKeep]
            if (exists("toKeepSR",inherits = FALSE)) { # must update toKeepSR in case that is a second reason to stop event
              toKeepSR <- toKeepSR[max_cellsKeep]
            }
            rm(whichID)
          }
          
        } else {
          if (size >= max_cells) {
            events <- NULL
          }
        }
        
      } else {
        events <- NULL
      }
      
      # drop or keep fire_start_location
      if (is.na(continue_to_burn_prob) | continue_to_burn_prob == 0L) {
        fire_start_location <- NULL
      } else {
        if (in.range(continue_to_burn_prob)) {
          fire_start_location <- fire_start_location[runif(length(fire_start_location)) <= continue_to_burn_prob]
        } else {
          # here is were we would handle methods for raster* or functions
          stop("Unsupported type: continue_to_burn_prob")
        }
      }
      
    fire_start_location <- c(fire_start_location, events)
    
    }
    
    # Convert the data back to raster
    wh <- spreads > 0
        if (return_table) {
          completed <- which(wh) %>%
            data.table::data.table(indices = ., id = spreads[.], active = FALSE)
          if (NROW(potentials) > 0) {
            active <- data.table::data.table(indices = potentials[, 2L],
                                 id = spreads[potentials[, 1L]],
                                 active = TRUE)
          } else {
            active <- data.table::data.table(indices = numeric(0), id = numeric(0),
                                 active = logical(0))
          }
      }
  
    
    
    if (return_table) {
        allCells <- data.table::rbindlist(list(completed, active))
        initEventID <- allCells[indices %in% initialfire_start_location, id]
        if (!all(is.na(initialfire_start_location))) {
          dtToJoin <- data.table(id = sort(initEventID), initialLocus = initialfire_start_location)
        } else {
          dtToJoin <- data.table(id = numeric(0), initialLocus = numeric(0))
        }
        if (fire_historyExists) {
          fire_historyInitialfire_start_location <- fire_history[, list(id = unique(id),
                                                       initialLocus = unique(initialLocus))]
          dtToJoin <- rbindlist(list(fire_historyInitialfire_start_location, dtToJoin))
        }
        data.table::setkeyv(dtToJoin, "id")
        data.table::setkeyv(allCells, "id")
        
        allCells <- dtToJoin[allCells]
        allCells[]
        return(allCells)
    }
    
    spre <- raster(habitat)
    spre[] <- 0
    spre[wh] <- spreads[wh]
      if (exists("potentials"))
        if (NROW(potentials) > 0)
          spre[potentials[, 1L]] <- spreads[potentials[, 2L]]
    
  return(spre)
}


in.range <-function (x, a = 0, b = 1) 
{
  if (is.null(x)) 
    return(NULL)
  if (!is.numeric(x)) {
    if (is(x, "Raster")) {
      x <- getValues(x)
    }
    else {
      stop("x must be numeric.")
    }
  }
  if (!is.numeric(a) || !is.numeric(b)) 
    stop("invalid (non-numeric) bounds.")
  if (is.na(a) || is.na(b)) 
    stop("invalid (NA) bounds.")
  if (a >= b) 
    stop("a cannot be greater than b.")
  return((x - a) * (b - x) >= 0)
}