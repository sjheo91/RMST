#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector RMST_ENS_cpp(DataFrame data, String method, double tau, NumericVector weight, double tol) {

  Function RMST_cpp( "RMST_cpp" ) ; 
  
  NumericVector rmst = rep(0.0, 3);
  NumericVector rmst0;
  NumericVector rmst1;
  NumericVector rmst2;
  
  if(method=="MCG"){
    rmst0 = RMST_cpp(data, method, tau, 2/3, 3, tol);
    rmst1 = RMST_cpp(data, method, tau, 2, 3, tol);
    rmst2 = RMST_cpp(data, method, tau, 6, 3, tol);
  }else{
    stop("Check method");
  }
  
  double rmst_ens = sum(rmst*weight)/sum(weight);
  return rmst_ens;
}
