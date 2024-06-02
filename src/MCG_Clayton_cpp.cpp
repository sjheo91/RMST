#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame MCG_Clayton_cpp(NumericVector time, NumericVector status, double theta) {
  
  if(theta<1e-09) { theta = 1e-09; }
  NumericVector u_t = unique(time).sort();
  int n = u_t.length();
  
  NumericVector m = Rcpp::rep(0.0, n);
  NumericVector d = Rcpp::rep(0.0, n);
  NumericVector c = Rcpp::rep(0.0, n);
  NumericVector R = Rcpp::rep(0.0, n);
  NumericVector A = Rcpp::rep(0.0, n);
  NumericVector S = Rcpp::rep(0.0, n);
  
  for(int i = 0; i < n; ++i){
    R[i] = time.length()-sum(m);
    m[i] = sum(ifelse(time==u_t[i], 1.0, 0.0));
    d[i] = sum(status*ifelse(time==u_t[i], 1.0, 0.0));
    c[i] = sum(ifelse(time==u_t[i], 1.0, 0.0)) - sum(status*ifelse(time==u_t[i], 1.0, 0.0));
    A[i] = pow((R[i]-d[i])/time.length(),-theta) - pow(R[i]/time.length(),-theta);
    if(abs(A[i])==1.0/0.0) { A[i] = 0.0; } 
    S[i] = pow(1 + sum(A),-1/theta);
  }
  
  DataFrame output = DataFrame::create(Named("time") = u_t, Named("surv") = S);
  return output;
}
