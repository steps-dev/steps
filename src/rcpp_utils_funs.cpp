#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<int> fast_match(IntegerVector x, IntegerVector y) {
  LogicalVector ind = in(x, y);
  int n = ind.size();
  std::vector<int> output;
  output.reserve(n);
  for (int i=0; i < n; ++i) {
    if (ind[i]) output.push_back(i+1);
  }
  return output;
}

// [[Rcpp::export]]
NumericVector pmax_zero(NumericVector X) {
  NumericVector Y = pmax(0, X);
  return Y;
}

// // [[Rcpp::export]]
// NumericVector rtnorm(NumericVector mu, NumericVector sd, double min, double max) {
//   // Obtain vector sizes
//   int n_mu = mu.size();
//   int n_sd = sd.size();
//   
//   // Check both vectors have elements
//   if(n_mu <= 0 || n_sd <= 0) {
//     Rcpp::stop("Both `mu` and `sd` inputs must have at least 1 element.");
//   }
//   
//   // Compare the vectors and recycle
//   NumericVector mu_vec = mu;
//   NumericVector sd_vec = sd;
//   if(n_mu > n_sd) {
//     mu_vec = mu;
//     sd_vec = rep_len(sd, n_mu);
//   } else if (n_sd > n_mu) {
//     mu_vec = rep_len(mu, n_sd);
//     sd_vec = sd;
//   }
//   
//   // Set n equal to the maximum length of mu or sd
//   int n = std::max(n_mu, n_sd);
//   
//   // Initiate output variable
//   NumericVector out (n);
//   
//   // Draw values with rejection if out of specified range
//   for (int i=0; i < n; ++i) {
//     bool valid = false;
//     double mu_draw = mu_vec[i];
//     double sd_draw = sd_vec[i];
//     while (!valid) {
//       NumericVector draw = Rcpp::rnorm(1, mu_draw, sd_draw);
//       if ((draw[0] <= max) & (draw[0] >= min)) {
//         out[i] = draw[0];
//         valid = true;
//       }
//     }
//   }
//   return(out) ;
// }

  