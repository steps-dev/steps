#' Create a growth transition function
#'
#' @description A growth transition function defines how spatial objects influence survival and fecundity. 
#' Functions that use elements of the landscape object to modify survival and fecundity are built-in for the user to select, however, a user may also provide a custom written growth function.
#' 
#' @rdname transition_function
#'
#' @param transition_matrix a symmetrical age-based (Leslie) or stage-based
#'   population structure matrix.
#' @param survival_layer a spatial object within the landscape object used to modify survival values.
#' @param fecundity_layer a spatial object within the landscape object used to modify fecundity values.
#' 
#' @return An object of class \code{transition_function}
#' 
#' @export
#'
#' @examples
#' 
#' test_transition_function <- modified_transition(egk_mat)

modified_transition <- function(transition_matrix,
                                survival_layer = "suitability",
                                fecundity_layer = "suitability") {
  
    idx <- which(transition_matrix != 0)
    is_recruitment <- upper.tri(transition_matrix)[idx]
    
    surv_vals <- transition_matrix[idx[!is_recruitment]]
    fec_vals <- transition_matrix[idx[is_recruitment]]
    
    dim <- nrow(transition_matrix)
  
  fun <- function (landscape, timestep) {
    
    # pull out survival/fecundity multipliers
    cell_idx <- which(!is.na(raster::getValues(landscape$population[[1]])))
    
    if (raster::nlayers(landscape$suitability) > 1) {
      surv_mult <- landscape[[survival_layer]][[timestep]][cell_idx]
      fec_mult <- landscape[[fecundity_layer]][[timestep]][cell_idx]
    } else {
      surv_mult <- landscape[[survival_layer]][cell_idx]
      fec_mult <- landscape[[fecundity_layer]][cell_idx]
    }
    
    # get per-cell versions of fecundity and survival values
    survs <- kronecker(surv_mult, t(surv_vals), "*")
    fecs <- kronecker(fec_mult, t(fec_vals), "*")
    
    # empty transition array to fill
    n_cells <- length(cell_idx)
    transition_array <- array(0, c(dim, dim, n_cells))
    
    # convert index from matrix to array
    addition <- dim ^ 2 * (seq_len(n_cells) - 1)
    idx_full <- as.numeric(outer(idx, addition, FUN = "+"))
    
    # put the surv/fec values back in (is_recruitment is recycled to match length)
    transition_array[idx_full[!is_recruitment]] <- survs
    transition_array[idx_full[is_recruitment]] <- fecs
    
    transition_array
    
  }

  as.transition_function(fun, info = list(transition_matrix = transition_matrix,
                                          survival_layer =  survival_layer,
                                          fecundity_layer = fecundity_layer))
  
}

# #' @rdname transition_function
# #'
# #' @param x an object to print or test as a transition_function object
# #' @param ... further arguments passed to or from other methods
# #'
# #' @export
# #'
# #' @examples
# #'
# #' print(test_transition_function)
# 
# print.transition_function <- function (x, ...) {
#   cat("This is a transition_function object")
# }

##########################
### internal functions ###
##########################

as.transition_function <- function (transition_function, info = NULL) {
  as_class(transition_function, "transition_function", "function", info = info)
}

print.dispersal_function <- function (x, ...) {
  print_info(x)
}
