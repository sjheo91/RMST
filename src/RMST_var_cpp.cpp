#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double RMST_var_cpp(DataFrame data, String method, double tau, double theta=1, int family=3, int n_boot=999, bool ensemble=false, NumericVector weight=NumericVector::create()) {
  
  Function subset( "[.data.frame" );
  Function RMST_cpp( "RMST_cpp" );
  Function RMST_ENS_cpp( "RMST_ENS_cpp" );
  
  int n = data.nrow();
  IntegerVector x = seq(1, n);
  NumericVector rmst;
  NumericVector rmst_vec = rep(0.0, n_boot);
  
  NumericVector v = NumericVector(n_boot*n, NumericVector::get_na()); 
  NumericMatrix ind = NumericMatrix(n_boot, n, v.begin());
  
  for(int i = 0; i < n_boot; ++i){
    ind(i,_) = sample(as<NumericVector>(x), n, true);
  }
  
  for(int i = 0; i < n_boot; ++i){
    if(ensemble==false){
      rmst = RMST_cpp(subset(data, ind(i,_), R_MissingArg), method, tau, theta, family);
    }else{
      rmst = RMST_ENS_cpp(subset(data, ind(i,_), R_MissingArg), method, tau, weight);
    }
    rmst_vec[i] = rmst[0];
  }
  
  return var(rmst_vec);
}
