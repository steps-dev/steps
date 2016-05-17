#include "RcppArmadillo.h"
using namespace Rcpp;

//'C++ connectivity function
//' @param dist Distances between patches (symetrical matrix)
//' @param alpha Exponential decay rate of patch connectivity (dispersion parameter)
//' @param beta double parameter that represents the shape of the dispersal kernel.
//' @param hanski_dispersal_kernal char if 'H' uses hanski(1994), if 'S' uses shaw(1995).
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix meta_dispersal_fun(NumericMatrix dist, double alpha, 
                                 double  beta=1, char hanski_dispersal_kernal = 'H') {
  arma::mat dist1 = as<arma::mat>(dist);
  arma::mat disp_mat(dist1.n_rows,dist1.n_cols);
  if(hanski_dispersal_kernal == 'H') disp_mat = exp(-alpha * dist1);
  if(hanski_dispersal_kernal == 'S') disp_mat = 1/(1+(alpha * arma::pow(dist1,beta)));
  return(wrap(disp_mat));
}

//'C++ colonisation function
//' @param s NumericVector that represents connectivity
//' @param y double colonisation parameter
//' @param c double colonisation scale parameter see Ovaskainen 2002.
//' @param col_meth char which method to use? 'H' = Hanski 1994, 'M'= Moilanen 2004, 'O'= Ovaskainen 2002. See decription for details.
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector meta_colonisation_fun(NumericVector s, double y, double c=1, char col_meth='H') {
  NumericVector cl(s.size()); 
  if(col_meth == 'H') cl = Rcpp::pow(s,2)/(Rcpp::pow(s,2) + (y*y));
  if(col_meth == 'M') cl = 1-Rcpp::exp(-y*s);
  if(col_meth == 'O') cl = Rcpp::pow(s,2)/(Rcpp::pow(s,2) + (1/c));
  return(cl);
}

// C++ Extinction function

//'C++ metapopulation function for a single timestep.
//' @param presence NumericVector Initial occupancy of each patch
//' @param dist_mat Exponential decay rate of patch connectivity (dispersion parameter)
//' @param Ei patch extinction rate at time i. note: In the future I need to pull this from demographic model.
//' @param y double colonisation parameter
//' @param c double colonisation scale parameter see Ovaskainen 2002.
//' @param col_meth char which method to use? 'H' = Hanski 1994, 'M'= Moilanen 2004, 'O'= Ovaskainen 2002. See decription for details.
//' @export
// [[Rcpp::export]]
NumericVector metapop(NumericVector presence, NumericMatrix dist_mat, NumericVector Ei, 
                      double y, double c = 1, char col_meth='H'){
  int presences = presence.length();
  arma::vec pres = presence;
  NumericVector s(presences);
  arma::mat dist_mat1 = as<arma::mat>(dist_mat);
  IntegerVector id = ifelse(presence>0,1,0);
  arma::vec ids = as<arma::vec>(id);
  arma::mat dm = dist_mat1.cols(arma::find(ids==1));
  NumericMatrix dm1 =wrap(dm);
  for(int i=0; i < presences; i++) {
    s[i] = sum(dm1(i,_));
  }
  NumericVector pa(presences); 
  NumericVector cl =  meta_colonisation_fun(s,y,c,col_meth);
  for (int j=0; j < presences; j++) {
    if (presence[j] == 0 && R::runif(0,1)  < cl[j])
      presence[j] = 1;
    else if (presence[j] == 1 && R::runif(0,1) < ((1 - cl[j]) * Ei[j]))
      presence[j] = 0;
  }
  return(presence);
}

//'Simulate a metapopulation system in C++
//'this function is not complete.
//' @param time Number of time steps
//' @param dist Distances between patches (symetrical matrix)
//' @param area Area of patches - This needs to be calculated somehow - using occupancy models?
//' @param presence Initial occupancies of patches. Must be presence 1 or absence 0.
//' @param y incidence function parameters
//' @param x incidence function parameters
//' @param e Minimum area of patches
//' @param alpha Exponential decay rate of patch connectivity (dispersion parameter)
//' @param beta double parameter that represents the shape of the dispersal kernel.
//' @param hanski_dispersal_kernal char if 'H' uses hanski(1994), if 'S' uses shaw(1995).
//' @param locations NULL or NumericMatrix Longitudes and latitudes of coordinates of the patches
//' @param c double colonisation scale parameter see Ovaskainen 2002.
//' @param col_meth char which method to use? 'H' = Hanski 1994, 'M'= Moilanen 2004, 'O'= Ovaskainen 2002. See decription for details.
//' @export
// [[Rcpp::export]]
NumericMatrix metapop_n(int time, NumericMatrix dist, NumericVector area, NumericVector presence, 
                 double y = 1, double x = 1, double e=1, double alpha = 1, double beta = 1, char hanski_dispersal_kernal = 'H',
                 Rcpp::Nullable<Rcpp::NumericMatrix> locations = R_NilValue, double c = 1, char col_meth='H'){
 arma::mat dist_mat = as<arma::mat>(meta_dispersal_fun(dist,alpha,beta,hanski_dispersal_kernal));
 dist_mat.diag().zeros();
 int ps = presence.size();
 arma::mat dist_mat2;
 dist_mat2.copy_size(dist_mat);
 for (int i=0; i < ps; i++) {
   dist_mat2.col(i) = dist_mat.col(i) * area[i];
 }
 NumericMatrix dist_mat3 = wrap(dist_mat2);
 int presences = presence.length();
 NumericMatrix presence_mat(presences,time+1);
 NumericVector E = e/Rcpp::pow(area,x);
 NumericVector Ei = ifelse(E > 1, 1, E);
 presence_mat(_,0) = presence;
 for (int i=0; i<time;i++){
   NumericVector tmp = presence_mat(_,i);
   NumericVector tmp1 = metapop(tmp,dist_mat3, Ei, y, c, col_meth);
   presence_mat(_,i+1)=tmp1;
 }
 return(presence_mat);
}

//'Simulate a metapopulation system in C++
//'this function is not complete.
//' @param nrep Number of simulations.
//' @param time Number of time steps
//' @param dist Distances between patches (symetrical matrix)
//' @param area Area of patches - This needs to be calculated somehow - using occupancy models?
//' @param presence Initial occupancies of patches. Must be presence 1 or absence 0.
//' @param y incidence function parameters
//' @param x incidence function parameters
//' @param e Minimum area of patches
//' @param alpha Exponential decay rate of patch connectivity (dispersion parameter)
//' @param beta double parameter that represents the shape of the dispersal kernel.
//' @param hanski_dispersal_kernal char if 'H' uses hanski(1994), if 'S' uses shaw(1995).
//' @param locations NULL or NumericMatrix Longitudes and latitudes of coordinates of the patches
//' @param c double colonisation scale parameter see Ovaskainen 2002.
//' @param col_meth char which method to use? 'H' = Hanski 1994, 'M'= Moilanen 2004, 'O'= Ovaskainen 2002. See decription for details.
//' @export
// [[Rcpp::export]]
List metapop_n_cpp(int nrep, int time, NumericMatrix dist, NumericVector area, NumericVector presence,
                   double y = 1, double x = 1, double e=1, double alpha = 1, double beta = 1, char hanski_dispersal_kernal = 'H',
                   Rcpp::Nullable<Rcpp::NumericMatrix> locations = R_NilValue, double c = 1, char col_meth='H'){
  List TempList(nrep);
   for (int ii=0;ii<nrep;ii++) {
      TempList[ii] = metapop_n(time = time, dist = dist, area = area, presence = presence,
                               y = y, x = x, e=e, alpha = alpha, beta = beta, 
                               hanski_dispersal_kernal = hanski_dispersal_kernal, locations = locations, c=c,col_meth=col_meth);
    }
    return(TempList);
}
