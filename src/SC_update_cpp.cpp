#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector SC_update_cpp(NumericVector time, NumericVector status, NumericVector Sx, NumericVector Sy, int family, double theta, double tol) {
  
  int n = time.length();
  NumericVector Sx_old = Sx * 1.0;
  NumericVector Sx_new = Sx * 1.0;
  NumericVector diff = Rcpp::rep(1.0, n);
  double RMSE = 1.0;
  
  while(RMSE > tol){
    for(int i = 0; i < n; ++i){
      NumericVector copu_cens = pow(pow(Sx_old[i],-theta)+pow(Sy,-theta)-1,-1/theta-1)/
        pow(pow(Sx_old,-theta)+pow(Sy,-theta)-1,-1/theta-1);
      copu_cens = ifelse(copu_cens==1.0/0.0, 0, copu_cens);
      NumericVector at_risk = ifelse(time>time[i], 1.0, 0.0);
      NumericVector non_at_risk = ifelse(time<=time[i], 1.0, 0.0);
      double S = (sum(at_risk)+sum(non_at_risk*(1.0-status)*copu_cens))/n;
      if(std::isnan(S)) S = 0.0;
      Sx_new[i] = S;
    }
    diff = Sx_old-Sx_new;
    RMSE = mean(pow(diff,2));
    Sx_old = Sx_new * 1.0;
  }
  return Sx_old;
}
