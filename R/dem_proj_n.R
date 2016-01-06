#' Demographic projection
#' 
#' @param v0 
#' @param tmat
#' @param matsd
#' @param stmat
#' @param estamb
#' @param estdem
#' @param equalsign
#' @return list of demograph projections
#' @export

demo_proj_n <- function (v0, tmat, matsd = NULL, stmat = NULL,
                         estamb = FALSE, estdem = FALSE, 
                         equalsign = TRUE, fecundity1 = FALSE, nrep = 10, 
                         time = 10)#,
          #management = NULL, round = TRUE) 
{
  vn <- NULL
  # vm <- NULL
  for (i in 1:nrep) {
    vn[[i]] <- cbind(v0, v0)
  }
  # call c++ function that does this loop.
  v <- demo_proj_n_cpp(vn, tmat, matsd = matsd, estamb = estamb, estdem = estdem, 
                       equalsign = equalsign, stmat = stmat, fecundity1 = fecundity1, nrep = nrep, 
                       time = time)
  vn <- lapply(v, function(x) x[, -1])
  vn <- list(vn = vn, mat = mat)
  return(vn)
}
