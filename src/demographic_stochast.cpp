#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
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