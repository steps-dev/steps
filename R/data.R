#' Eastern Grey Kangaroo habitat suitability
#'
#' A raster layer containing the predicted relative habitat suitability
#' for the Eastern Grey Kangaroo.
#'
#' @format A raster layer with 35 rows, 36 columns, and 1260 cells:
#' \describe{
#'   \item{resolution}{size of cells in meters}
#'   \item{extent}{cartesian boundaries in meters}
#'   \item{coord. ref.}{spatial projection of layer}
#' }
"egk_hab"

#' Eastern Grey Kangaroo populations
#'
#' A raster stack containing intital populations for each life-stage
#' of the Eastern Grey Kangaroo.
#'
#' @format A raster stack with 35 rows, 36 columns, 1260 cells, and 3 layers:
#' \describe{
#'   \item{resolution}{size of cells in meters}
#'   \item{extent}{cartesian boundaries in meters}
#'   \item{coord. ref.}{spatial projection of layer}
#' }
"egk_pop"

#' Eastern Grey Kangaroo carrying capacity
#'
#' A raster stack containing the total number of Eastern Grey Kangaroos
#' each grid cell can support.
#'
#' @format A raster layer with 35 rows, 36 columns, and 1260 cells:
#' \describe{
#'   \item{resolution}{size of cells in meters}
#'   \item{extent}{cartesian boundaries in meters}
#'   \item{coord. ref.}{spatial projection of layer}
#' }
"egk_k"

#' Eastern Grey Kangaroo life-stage transition matrix
#'
#' A matrix containing the survival and fecundity of Eastern Grey Kangaroos
#' at each of three life-stages - juvenile, subadult, and adult.
#'
#' @format A matrix with 3 rows and 3 columns:
#' \describe{
#'   \item{juvenile}{this column specifies the survival of juveniles into the next stage}
#'   \item{subadult}{this column specifies the survival of subadults into the next stage}
#'   \item{adult}{this column specifies the survival of adults into the next year and the number of offsring produced (first row)}
#' }
"egk_mat"

#' Eastern Grey Kangaroo state object
#'
#' An object containing the habitat, demographic and populations of Eastern Grey Kangaroos.
#' This is created with the \link[steps]{state} function.
#'
#' @format A state object with:
#' \describe{
#'   \item{population}{spatial layers describing the population}
#'   \item{habitat}{spatial layers describing habitat suitability and carrying capacity}
#'   \item{demography}{matrices describing survival and fecundity}
#' }
"egk_state"

#' Eastern Grey Kangaroo habitat disturbance
#'
#' A raster stack containing values for modifying the habitat.
#'
#' @format A raster stack with 35 rows, 36 columns, 1260 cells, and 20 layers:
#' \describe{
#'   \item{resolution}{size of cells in meters}
#'   \item{extent}{cartesian boundaries in meters}
#'   \item{coord. ref.}{spatial projection of layer}
#' }
"egk_dist"