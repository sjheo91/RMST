#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame MCG_cpp(NumericVector time, NumericVector status, double theta, int family) {
  
  if(family==3){
    DataFrame output = MCG_Clayton_cpp(time, status, theta);
    return output;
  }else if(family==4){
    DataFrame output = MCG_Gumbel_cpp(time, status, theta);
    return output;
  }else if(family==5){
    DataFrame output = MCG_Frank_cpp(time, status, theta);
    return output;
  }else{
    stop("Check family");
  }
}
