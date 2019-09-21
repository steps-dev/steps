#' Eastern Grey Kangaroo example data
#'
#' Example data for simulating spatial population dynamics of Eastern Grey
#' Kangaroos in a hypothetical landscape.
#' 
#' \describe{
#'   \item{egk_hab}{A raster layer containing the predicted relative habitat
#'   suitability for the Eastern Grey Kangaroo.}
#'   \item{egk_pop}{A raster stack containing initial populations for each
#'   life-stage of the Eastern Grey Kangaroo.}
#'   \item{egk_k}{A raster layer containing the total number of Eastern Grey
#'   Kangaroos each grid cell can support.}
#'   \item{egk_mat}{A matrix containing the survival and fecundity of Eastern
#'   Grey Kangaroos at each of three life-stages - juvenile, subadult, and adult.}
#'   \item{egk_mat_stoch}{A matrix containing the uncertainty around survival
#'   and fecundity of Eastern Grey Kangaroos at each of three life-stages -
#'   juvenile, subadult, and adult.}
#'   \item{egk_sf}{A raster stack containing values for modifying survival and
#'   fecundities - each is raster is named according to the timestep and position
#'   of the life-stage matrix to be modified.}
#'   \item{egk_fire}{A raster stack containing values for modifying the habitat
#'   - in this case the proportion of landscape remaining after fire.}
#'   \item{egk_origins}{A raster stack containing locations and counts of where
#'   to move individual kangaroos from.}
#'   \item{egk_destinations}{A raster stack containing locations and counts of where to
#'   move individual kangaroos to.}
#'   \item{egk_road}{A raster stack containing values for modifying the habitat
#'   - in this case the proportion of habitat remaining after the construction of a road.}
#' }
#' @format Misc data
#' @name egk
"egk_hab"

#' @rdname egk
"egk_pop"

#' @rdname egk
"egk_k"

#' @rdname egk
"egk_mat"

#' @rdname egk
"egk_mat_stoch"

#' @rdname egk
"egk_sf"

#' @rdname egk
"egk_fire"

#' @rdname egk
"egk_origins"

#' @rdname egk
"egk_destinations"

#' @rdname egk
"egk_road"