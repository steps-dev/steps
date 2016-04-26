#' @useDynLib ramasr
#' @importFrom Rcpp sourceCpp
NULL

#' demo_proj_n
#' 
#' Demographic projection
#' 
#' Make deterministic and stochastic demographic projections according to a transition matrix.
#' 
#' @param v0 vector Starting abundances at each timestep.
#' @param tmat matrix. Transition matrix
#' @param matsd matrix. Matrix with the standard deviation of the probabilities in tmat. 
#' @param stmat martix. Matrix indicating for each transition probability in mat which part (i.e. which proportion) should be considered resulting from fecundity.
#' @param estamb Logical. Should environmental stochasticity be considered to projet the dynamics of the population?
#' @param estdem Logical. Should demographic stochasticity be employed to project the dynamics of the population?
#' @param equalsign Logical. Should the environmental deviations have all the same sign and magnitude? 
#' @param tmat_fecundity Logical. Should the first row of tmat as fecundities?
#' @param nrep int number of simulations to run.
#' @param time int length of the demographic trajectory.
#' @return list of demograph projections
#' @author Skipton Woolley
#' @seealso \code{link{demo_proj_n_cpp}}
#' @examples 
#' tmat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
#' v0 <- c(80,20,0)
#' sim_t10_rep100 <- demo_proj_n(v0=v0,tmat=tmat,time=10,nrep=100) 
#' @export
demo_proj_n <- function (v0, tmat, matsd = NULL, stmat = NULL,
                         estamb = FALSE, estdem = FALSE, 
                         equalsign = TRUE, tmat_fecundity = FALSE, nrep = 10, 
                         time = 10)
{
  vn <- NULL
  # vm <- NULL
  for (i in 1:nrep) {
    vn[[i]] <- base::cbind(v0, v0)
  }
  # call c++ function that does this loop.
  v <- ramasr::demo_proj_n_cpp(vn, tmat, matsd = matsd, estamb = estamb, estdem = estdem, 
                       equalsign = equalsign, stmat = stmat, tmat_fecundity = tmat_fecundity,
                       nrep = nrep, time = time)
  vn <- base::lapply(v, function(x) x[,-1])
  vn <- base::list(vn = vn, tmat = tmat)
  return(vn)
}
