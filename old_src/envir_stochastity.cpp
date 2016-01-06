// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::mat envir_stochast(arma::mat tmat, arma::mat sdmat, bool equalsign = true)
 {
   arma::vec tmat_v = arma::vectorise(tmat);
   arma::vec sdmat_v = arma::vectorise(sdmat);
   int nvals = tmat_v.size(); 
   arma::vec mat_v = tmat_v.zeros();
   arma::vec deriv_v = tmat_v.zeros();
   int nc = tmat.n_cols;
   int nr = tmat.n_rows;     
   arma::mat mat1; 
   if (equalsign == false) {
       for(int i = 0; i<nvals; i++){
          mat_v[i] = R::rnorm(tmat_v[i], sdmat_v[i]);
       }
       mat1.insert_cols(0, mat_v);
       mat1.reshape(nr, nc);
   }     
   if (equalsign == true) {
     for(int i = 0; i<nvals; i++){
            Rcpp::RNGScope tmp;
            double draw = R::runif(-1,1);
       if (draw >= 0) {
         draw = 1;
       } else {
         draw = -1;
       }
       // Rcpp::Rcout << R::rnorm(0, sdmat_v[i])*draw << std::endl;
       deriv_v[i] = R::rnorm(0, sdmat_v[i])*draw;
       }                     
           deriv_v = abs(deriv_v);
           // Rcpp::Rcout << deriv_v << std::endl;
           mat1.insert_cols(0, deriv_v);
           mat1.reshape(nr, nc);
           // Rcpp::Rcout << mat1 << std::endl;
           mat1 = tmat + mat1;
           
    }       
    // mat[mat < 0] = 0;
    // mat[-1, ][mat[-1, ] > 1] = 1;
    return(mat1);
}
   