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
NumericVector demographic_stochast(NumericVector v, NumericMatrix tmat)//, NumericMatrix stmat = NULL,  Logical fecundity1 = TRUE) 
{
  int tmncols = tmat.ncol();
  NumericMatrix sij(tmncols-1,tmncols);
  NumericVector sj(tmncols-1);
  NumericVector result(tmncols);
  int bs = sum(v);
  NumericVector Bi(bs);
    
    int b = 0;
    for(int i = 0; i<tmncols; i++) {
        if (tmat(0, i) > 0) 
          for(int ii = 0; ii<v[i]; ii++){
              Bi[ii+b] = R::rpois(tmat(0, i));
              //Rcpp::Rcout << ii << std::endl;
          }
        for (int j = 0; j<tmncols-1; j++) {
          if (tmat(j + 1, i) > 0) {
           if (tmat(j + 1, i) > 1) 
               sj[j] = sum(rpois(v[i], tmat(j + 1, i)));
               else sj[j] = sum(rbinom(v[i], 1, tmat(j + 1, i)));
            }
          else sj[j] = 0;
        }
      sij(_,i) = sj;
      b = b+v[i];
    }  
    result[0] = sum(Bi);
    for(int k = 0; k<tmncols-1;k++){
      result[k+1] = sum(sij(k,_));
    }
    return(result);
}