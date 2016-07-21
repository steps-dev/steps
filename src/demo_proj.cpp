 // we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;

//'Demographic stochastic function in C++
//' @param v NumericVector. Vector with the initial abundance of each stage.
//' @param tmat NumericMatrix. Transition matrix.
//' @param stmat NumericMatrix. Matrix indicating for each transition probability in mat which part (i.e. which proportion) should be considered resulting from fecundity.
//' @paran tmat_fecundity Logical. If TRUE use first row of transition matrix as fecunity. 
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector demographic_stochast(NumericVector v, 
                                   NumericMatrix tmat,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> stmat = R_NilValue,
                                   bool tmat_fecundity = true) 
{
  int tmncols = tmat.ncol();
  NumericMatrix sij(tmncols-1,tmncols);
  NumericVector sj(tmncols-1);
  NumericVector result(tmncols);
  int bs = sum(v);
  NumericVector Bi(bs);
  
  // if (stmat.isNull() && tmat_fecundity == true) {
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
  // }
  return(result);
}


//'environmental stochastic function in C++
//' @param tmat NumericMatrix. Vector with the initial abundance of each stage.
//' @param matsd NumericMatrix. Transition matrix.
//' @param equalsign bool. Should the environmental deviations have all the same sign and magnitude? TRUE or FALSE
//' @export
// [[Rcpp::export]]
NumericMatrix envir_stochast(NumericMatrix tmat, NumericMatrix matsd, bool equalsign = true)
{ 
  arma::mat tmat1 = as<arma::mat>(tmat);
  arma::mat matsd1 = as<arma::mat>(matsd);
  arma::vec tmat_v = arma::vectorise(tmat1);
  arma::vec matsd_v = arma::vectorise(matsd1);
  int nvals = tmat_v.size(); 
  arma::vec mat_v = tmat_v.zeros();
  arma::vec deriv_v = tmat_v.zeros();
  int nc = tmat1.n_cols;
  int nr = tmat1.n_rows;     
  arma::mat mat1; 
  if (equalsign == false) {
    for(int i = 0; i<nvals; i++){
      mat_v[i] = R::rnorm(tmat_v[i], matsd_v[i]);
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
      deriv_v[i] = R::rnorm(0, matsd_v[i])*draw;
    }                     
    deriv_v = abs(deriv_v);
    mat1.insert_cols(0, deriv_v);
    mat1.reshape(nr, nc);
    mat1 = tmat1 + mat1;
    }       
  return(wrap(mat1));
}

//' Single time step demographic projection function in C++
//' 
//' @param v0 NumericVector. Vector with the initial abundance of each stage.
//' @param tmat NumericMatrix. Transition matrix.
//' @param matsd NumericMatrix. Matrix with the standard deviation of the probabilities in tmat. 
//' @param stmat NumericMatrix. Matrix indicating for each transition probability in mat which part (i.e. which proportion) should be considered resulting from fecundity.
//' @param estamb bool. Environmental stochasticity included in population dynamics?
//' @param estdem bool. Demographic stochasticity included in population dynamics?
//' @param equalsign bool. Should the environmental deviations have all the same sign and magnitude?
//' @param tmat_fecundity bool. Should the first row of tmat as fecundities? TRUE or FALSE
//' @export
// [[Rcpp::export]]
NumericVector demo_proj(NumericVector v0, 
                        NumericMatrix tmat, 
                        Rcpp::Nullable<Rcpp::NumericMatrix> matsd= R_NilValue, 
                        Rcpp::Nullable<Rcpp::NumericMatrix> stmat = R_NilValue,
                        bool estamb=false,
                        bool estdem=false,
                        bool equalsign=true,
                        bool tmat_fecundity=false) 
  {
    NumericVector v1 = v0;
    arma::vec av = as<arma::vec>(v0);
    arma::mat m1;
    arma::mat em;
    m1.insert_cols(0, av);
    arma::mat tmat1 = as<arma::mat>(tmat);
    if (estamb == false && estdem == false)
      v1 = tmat1 * m1;
      wrap(v1);
      if (estamb == false && estdem == true && stmat.isNull() && tmat_fecundity == true) 
        v1 = demographic_stochast(v0, tmat);
        // if (estamb == false && estdem == true && stmat.isNull() && tmat_fecundity == false)
          // v1 =  demographic_stochast(v0, tmat, tmat_fecundity = false);
          // if (estamb == false && estdem == true && stmat.isNotNull()) 
          //   v1 =  demographic_stochast(v0, tmat, stmat = stmat);
            // if (estamb == true && estdem == false) {
            //   if (matsd.isNull()) 
            //     stop("there is not SD matrix provided\n (argument matsd=NULL)");
            //     em =  envir_stochast(tmat, matsd, equalsign = equalsign);
            //     v1 = v1 * m1;
            //     wrap(v1);
            // }
            // if (estamb == true & estdem == true) {
            //   if (matsd.isNull()) 
            //     stop("there is not SD matrix provided\n (argument matsd=NULL)");
            //     if (stmat.isNull() && tmat_fecundity == true)
            //       em = envir_stochast(tmat, matsd, equalsign = equalsign);
            //       v1 = demographic_stochast(v0, tmat = em);
            //       if (stmat.isNull() && tmat_fecundity == false)
            //         em = as<arma::mat>(envir_stochast(tmat, matsd, equalsign = equalsign));
            //         v1 = demographic_stochast(v0, tmat = wrap(em), tmat_fecundity = false);
            //         if (stmat.isNotNull()) 
            //           arma::mat em = as<arma::mat>(envir_stochast(tmat, matsd, equalsign = equalsign));
            //           v1 = demographic_stochast(v0, tmat = wrap(em), stmat = stmat);
            // }    
            return(v1);
}

//' Multiple time step and repetition demographic projection function in C++
//' 
//' @param vn List. List of vectors for simulation.
//' @param tmat NumericMatrix. Transition matrix.
//' @param matsd NumericMatrix. Transtion matrix error.
//' @param stmat NumericMatrix. Matrix indicating for each transition probability in mat which part (i.e. which proportion) should be considered resulting from fecundity
//' @param estamb bool. Environmental stochasticity included in population dynamics?
//' @param estdem bool. Demographic stochasticity included in population dynamics?
//' @param equalsign Logical. Should the environmental deviations have all the same sign and magnitude?
//' @param tmat_fecundity bool Should the first row of tmat as fecundities? TRUE or FALSE
//' @param nrep int number of simulations
//' @param time int number of time-steps.
//' @export
// [[Rcpp::export]]
List demo_proj_n_cpp(List vn, NumericMatrix tmat, Rcpp::Nullable<Rcpp::NumericMatrix> matsd= R_NilValue, 
                       Rcpp::Nullable<Rcpp::NumericMatrix> stmat = R_NilValue,
                       bool estamb=false, bool estdem=false, bool equalsign=true,  bool tmat_fecundity=false,
                       int nrep = 1, int time = 10)//, Rcpp::Nullable<Rcpp::NumericMatrix> management= R_NilValue, bool round = true) 
{
  List vn1 = vn;
  for (int i=0;i<time;i++) {
    for (int ii=0;ii<nrep;ii++) {
      NumericMatrix vii = vn[ii];
      NumericVector vii_i = vii(_,i+1);
      NumericVector v = demo_proj(vii_i, tmat, matsd, stmat, estamb, estdem,
                    equalsign, tmat_fecundity);
      arma::vec v1 = as<arma::vec>(v);
      arma::mat m1 = as<arma::mat>(vn1[ii]);
      m1.insert_cols(m1.n_cols, v1);
      vn1[ii] =  m1; 
      }
  }
  return(vn1);
}
