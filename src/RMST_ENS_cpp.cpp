#include <Rcpp.h>
using namespace Rcpp;

Function RMST_cpp( "RMST_cpp" ) ; 

// [[Rcpp::export]]
double RMST_ENS_cpp(DataFrame data, String method, double tau, NumericVector weight) {
  
  NumericVector rmst = rep(0.0, 3);
  double rmst0;
  double rmst1;
  double rmst2;
  
  if(method=="MCG"){
    rmst0 = RMST_cpp(data, method, tau, 2/3, 3);
    rmst1 = RMST_cpp(data, method, tau, 2, 3);
    rmst2 = RMST_cpp(data, method, tau, 6, 3);
  }else{
    stop("Check method");
  }
  
  double rmst_ens = sum(rmst*weight)/sum(weight);
  return rmst_ens;
}
