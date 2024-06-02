#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double RMST_ENS_cpp(DataFrame data, String method, double tau, NumericVector weight) {
  
  NumericVector rmst = rep(0.0, 3);
  
  if(method=="MCG"){
    rmst[0] = RMST_cpp(data, method, tau, 2/3, 3);
    rmst[1] = RMST_cpp(data, method, tau, 2, 3);
    rmst[2] = RMST_cpp(data, method, tau, 6, 3);
  }else{
    stop("Check method");
  }
  
  double rmst_ens = sum(rmst*weight)/sum(weight);
  return rmst_ens;
}
