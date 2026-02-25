#' Eastern Grey Kangaroo example data
#'
#' Example data for simulating spatial population dynamics of Eastern Grey
#' Kangaroos in a hypothetical landscape.
#' 
#' \describe{
#'   \item{egk_hab}{A function that generates the SpatRaster containing the
#'    predicted relative habitat
#'   suitability for the Eastern Grey Kangaroo.}
#'   \item{egk_pop}{A function that generates the SpatRaster stack containing
#'    initial populations for each
#'   life-stage of the Eastern Grey Kangaroo.}
#'   \item{egk_k}{A function that generates the SpatRaster layer containing the
#'    total number of Eastern Grey
#'   Kangaroos each grid cell can support.}
#'   \item{egk_mat}{A matrix containing the survival and fecundity of Eastern
#'   Grey Kangaroos at each of three life-stages - juvenile, subadult, and 
#'   adult.}
#'   \item{egk_mat_stoch}{A matrix containing the uncertainty around survival
#'   and fecundity of Eastern Grey Kangaroos at each of three life-stages -
#'   juvenile, subadult, and adult.}
#'   \item{egk_sf}{A function that generates the SpatRaster stack containing 
#'   values for modifying survival and
#'   fecundities - each is raster is named according to the timestep and 
#'   position
#'   of the life-stage matrix to be modified.}
#'   \item{egk_fire}{A function that generates the SpatRaster stack containing 
#'   values for modifying the habitat
#'   - in this case the proportion of landscape remaining after fire.}
#'   \item{egk_origins}{A function that generates the SpatRaster stack 
#'   containing locations and counts of where
#'   to move individual kangaroos from.}
#'   \item{egk_destinations}{A function that generates the SpatRaster stack
#'    containing locations and counts of where to
#'   move individual kangaroos to.}
#'   \item{egk_road}{A function that generates the SpatRaster stack containing
#'    values for modifying the habitat
#'   - in this case the proportion of habitat remaining after the construction
#'    of a road.}
#' }
#' @format Misc data
#' @name egk
#' @docType data
NULL

### Raster constructor functions ------------
#' @rdname egk
#' @export
egk_hab <- function() {
  terra::rast(system.file("extdata/egk_hab.tif", package = "steps"))
}

#' @rdname egk
#' @export
egk_pop <- function() {
  terra::rast(system.file("extdata/egk_pop.tif", package = "steps"))
}

#' @rdname egk
#' @export
egk_k <- function() {
  terra::rast(system.file("extdata/egk_k.tif", package = "steps"))
}

#' @rdname egk
#' @export
egk_sf <- function() {
  terra::rast(system.file("extdata/egk_sf.tif", package = "steps"))
}

#' @rdname egk
#' @export
egk_fire <- function() {
  terra::rast(system.file("extdata/egk_fire.tif", package = "steps"))
}

#' @rdname egk
#' @export
egk_origins <- function() {
  terra::rast(system.file("extdata/egk_origins.tif", package = "steps"))
}

#' @rdname egk
#' @export
egk_destinations <- function() {
  terra::rast(system.file("extdata/egk_destinations.tif", package = "steps"))
}

#' @rdname egk
#' @export
egk_road <- function() {
  terra::rast(system.file("extdata/egk_road.tif", package = "steps"))
}


### Lazy-loaded matrices ======

#' @rdname egk
"egk_mat"

#' @rdname egk
"egk_mat_stoch"