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