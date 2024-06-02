#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double RMST_ENS_cpp(DataFrame data, String method, double tau, NumericVector weight, double tol = 1e-6) {

  Function RMST_cpp( "RMST_cpp" ) ; 

  NumericVector rmst;
  NumericVector rmst_1;
  NumericVector rmst_2;
  NumericVector rmst_3;
  
  if(method=="MCG"|method=="SC"){
    rmst_1 = RMST_cpp(data, method, tau, 2/3, 3, tol);
    rmst_2 = RMST_cpp(data, method, tau, 2, 3, tol);
    rmst_3 = RMST_cpp(data, method, tau, 6, 3, tol);
    rmst = {rmst_1[0], rmst_2[0], rmst_3[0]};
  }else{
    stop("Check method");
  }
  
  double rmst_ens = sum(rmst*weight)/sum(weight);
  return rmst_ens;
}
