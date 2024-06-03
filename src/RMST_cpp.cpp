#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double RMST_cpp(DataFrame data, String method, double tau, int family, double theta, double tol = 1e-6) {

  Function MCG_cpp( "MCG_cpp" ) ; 
  Function SC_copula_cpp( "SC_copula_cpp" ) ; 
  Function sort_by( "sort_by" ) ; 
  
  double rmst;
  DataFrame ft;
  
  NumericVector time = sort_by(data["time"], data["time"]);
  NumericVector status = sort_by(data["status"], data["time"]);
  
  if(method=="MCG"){
    ft = MCG_cpp(time, status, family, theta);
  }else if(method=="SC"){
    ft = SC_copula_cpp(time, status, family, theta, tol);
  }else{
    stop("Check method");
  }
  
  NumericVector t = ft["time"];
  NumericVector s = ft["surv"];
  
  NumericVector t_sub = t[t<tau];
  NumericVector s_sub = s[t<tau];
    
  t_sub.insert(t_sub.begin(), 0);
  t_sub.insert(t_sub.end(), tau);
  s_sub.insert(s_sub.begin(), 1);
  
  rmst = sum(diff(t_sub)*s_sub);
  return rmst;
}
