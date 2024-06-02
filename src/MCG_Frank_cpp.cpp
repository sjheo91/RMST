// [[Rcpp::export]]
DataFrame MCG_Frank_cpp(NumericVector time, NumericVector status, double theta) {
  
  if(abs(theta)<1e-06) { theta = 1e-06; }
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
    A[i] = log((exp(-theta * (R[i] - d[i])/time.length()) - 1.0) / (exp(-theta * R[i]/time.length()) - 1.0));
    if(abs(A[i])==1.0/0.0) { A[i] = 0.0; } 
    S[i] = -1/theta * log(1 + (exp(-theta) - 1) * exp(sum(A)));
  }
  
  DataFrame output = DataFrame::create(Named("time") = u_t, Named("surv") = S);
  return output;
}
