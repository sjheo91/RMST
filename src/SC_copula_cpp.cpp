#include <Rcpp.h>
using namespace Rcpp;

Function SC_update_cpp( "SC_update_cpp" ) ; 

// [[Rcpp::export]]
DataFrame SC_copula_cpp(NumericVector time, NumericVector status, int family, double theta, double tol) {
  
  int n = time.length();
  
  NumericVector Se = Rcpp::rep(0.1, n);
  NumericVector Se_new = Rcpp::rep(0.1, n);
  
  NumericVector Sc = Rcpp::rep(0.1, n);
  NumericVector Sc_new = Rcpp::rep(0.1, n);
  
  NumericVector diff = Rcpp::rep(1.0, 2);
  double RMSE = 1.0;
  
  while(RMSE > tol){
    Se_new = SC_update_cpp(time, status, Se, Sc, family, theta, tol);
    Sc_new = SC_update_cpp(time, 1-status, Sc, Se_new, family, theta, tol);
    
    diff = (Se-Se_new);
    RMSE = mean(pow(diff,2));
    
    Se = Se_new*1.0;
    Sc = Sc_new*1.0;
  }
  
  DataFrame output = DataFrame::create(Named("time") = time, Named("surv") = Se);
  
  return output;
}