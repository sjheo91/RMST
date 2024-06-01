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

// [[Rcpp::export]]
DataFrame MCG_Gumbel_cpp(NumericVector time, NumericVector status, double theta) {
  
  if(theta<0) { theta = 0; }
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
    A[i] = pow(-log((R[i]-d[i])/time.length()),theta) - pow(-log(R[i]/time.length()),theta);
    if(abs(A[i])==1.0/0.0) { A[i] = 0.0; } 
    S[i] = exp(-pow(sum(A),1/theta));
  }
  
  DataFrame output = DataFrame::create(Named("time") = u_t, Named("surv") = S);
  return output;
}

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

// [[Rcpp::export]]
NumericVector sort_by(NumericVector x, NumericVector y) {
  IntegerVector idx = seq_along(x) - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
  return x[idx];
}

// [[Rcpp::export]]
NumericVector SC_update_cpp(NumericVector time, NumericVector status, 
                            NumericVector Sx, NumericVector Sy, int family, double theta, double tol) {
  
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

// [[Rcpp::export]]
DataFrame SC_copula_cpp(NumericVector time, NumericVector status, double theta, int family, double tol=1e-6) {
  
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

// [[Rcpp::export]]
double RMST_cpp(DataFrame data, String method, double tau, double theta, int family) {
  
  double rmst;
  DataFrame ft;
  
  NumericVector time = sort_by(data["time"], data["time"]);
  NumericVector status = sort_by(data["status"], data["time"]);
  
  if(method=="MCG"){
    ft = MCG_cpp(time, status, theta, family);
  }else if(method=="SC"){
    ft = SC_copula_cpp(time, status, theta, family);
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

// [[Rcpp::export]]
double RMST_var_cpp(DataFrame data, String method, double tau, double theta=1, int family=3, int n_boot=999, bool ensemble=false, NumericVector weight=NumericVector::create()) {
  
  Function subset( "[.data.frame" );
  
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
