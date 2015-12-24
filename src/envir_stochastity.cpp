// // -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// 
// // we only include RcppArmadillo.h which pulls Rcpp.h in for us
// #include "RcppArmadillo.h"
// using namespace Rcpp;
// // via the depends attribute we tell Rcpp to create hooks for
// // RcppArmadillo so that the build process will know what to do
// //
// // [[Rcpp::depends(RcppArmadillo)]]
// 
// // simple example of creating two matrices and
// // returning the result of an operatioon on them
// //
// // via the exports attribute we tell Rcpp to make this function
// // available from R
// //
// // [[Rcpp::export]]
// NumericMatrix envir_stochast(NumericMatrix tmat, NumericMatrix sdmat, bool equalsign = true)
// {
//   NumericMatrix mat = tmat;
//   arma::vec tmat_v = arma::vectorise(tmat);
//   arma::vec sdmat_v = arma::vectorise(sdmat);
//   int nvals = tmat_v.size(); 
//     if (equalsign == false) {
//     mat = R::rnorm(nvals, mean = tmat_v,sd = sdmat_v);
//   }
// //   if (equalsign == true) {
// //     deviates <- abs(rnorm(length(mat), mean = 0, sd = as.vector(matsd))) * sample(c(-1, 1), 1)
// //     mat <- mat + matrix(deviates, nrow = dim(mat)[1], ncol = dim(mat)[2])
// //   }
//   mat[mat < 0] = 0;
//   // mat[-1, ][mat[-1, ] > 1] = 1;
//   return(mat);
// }
//   