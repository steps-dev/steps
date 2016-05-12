#include "RcppArmadillo.h"
using namespace Rcpp;

//'C++ dispersal function
//' @param dist Distances between patches (symetrical matrix)
//' @param alpha Exponential decay rate of patch connectivity (dispersion parameter)
//' @param beta double parameter that represents the shape of the dispersal kernel.
//' @paran hanski_dispersal_kernal bool if true uses hanski(1994), if false uses shaw(1995).
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix meta_dispersal_fun(NumericMatrix dist, double alpha, 
                                 double  beta=1, bool hanski_dispersal_kernal = true) {
  arma::mat dist1 = as<arma::mat>(dist);
  arma::mat disp_mat(dist1.n_rows,dist1.n_cols);
  if(hanski_dispersal_kernal == true) disp_mat = exp(-alpha * dist1);
  if(hanski_dispersal_kernal == false) disp_mat = 1/(1+(alpha * arma::pow(dist1,beta)));
  return(wrap(disp_mat));
}

//'Simulate a metapopulation system in C++
//' @param time Number of time steps
//' @param dist Distances between patches (symetrical matrix)
//' @param area Area of patches - This needs to be calculated somehow - using occupancy models?
//' @param presence Initial occupancies of patches. Must be presence 1 or absence 0.
//' @param y incidence function parameters
//' @param x incidence function parameters
//' @param e Minimum area of patches
//' @param alpha Exponential decay rate of patch connectivity (dispersion parameter)
//' @param locations NULL or NumericMatrix Longitudes and latitudes of coordinates of the patches
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix sim_metapop(int time, NumericMatrix dist, NumericVector area, NumericVector presence, 
                 int y = 1, int x = 1, int e=1, double alpha = 1, double beta = 1, bool hanski_dispersal_kernal = true,
                 Rcpp::Nullable<Rcpp::NumericMatrix> locations = R_NilValue){
  arma::mat edis = as<arma::mat>(meta_dispersal_fun(dist,alpha,beta,hanski_dispersal_kernal));
  edis.diag().zeros(); 
  return(wrap(edis));
}
//   arma::mat pmat
//   edis <- sweep(edis, 2, A, "*")
//   E <- e/A^x
//   E <- ifelse(E > 1, 1, E)
//   if (locations.isNull())NumericVector locations = R::cmdscale(d);
//     pmat <- matrix(0, nrow = length(p), ncol = steps + 1)
//     pmat[, 1] <- p
//     for (i in 1:steps) pmat[, i + 1] <- metastep(pmat[, i], 
//          edis, E, y)
//       out <- list(p = pmat, d = edis, A = A, y = y, x = x, e = e, 
//                   alpha = alpha, locations = locations)
//       out$J.obs <- rowSums(pmat[, -1, drop = FALSE])/steps
//       out$P.obs <- colSums(pmat)
//       S <- rowSums(edis)
//       C <- S^2/(S^2 + y)
//       out$J.pot <- C/(C + E - C * E)
//       out$S.pot <- S
//       out$C.pot <- C
//       class(out) <- "metacycle"
//     out
// }
