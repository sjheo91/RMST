#' @export
#'
#' @import survival
#' @import dplyr
#' @importFrom foreach %dopar% foreach
#' @importFrom iterators icount
#' @importFrom doParallel registerDoParallel

RMST_var_cpp <- function(data, method, tau, theta, family, ensemble=False, tol=1e-6,
                     n_boots = 999,
                     n_cores = parallel::detectCores()) {
  
  # Construct cluster
  cl = parallel::makeCluster(n_cores)
  
  # After the function is run, close the cluster.
  on.exit(parallel::stopCluster(cl))
  
  # Register parallel backend
  doParallel::registerDoParallel(cl)
  
  # Compute estimates
  output <- foreach::foreach(i = iterators::icount(n_boots), # Perform n simulations
                                .combine = "rbind",
                                .packages = "Rcpp2doParallel") %dopar% {
    boot_data = data[sample(1:nrow(data), replace=T),]
    if(ensemble){
      result = RMST_cpp(boot_data, method, tau, theta, family, tol);
    }else{
      result = RMST_ENS_cpp(boot_data, method, tau, weight, tol);
    }
    result;
  }
  output
}

RMST_indep <- function(data, tau){
  
  ft <- survfit(Surv(data$time, data$status) ~ 1)
  
  ind <- which(ft$time<=tau)
  time <- ft$time[ind]
  surv <- ft$surv[ind]
  
  area <- diff(c(0,time,tau))*c(1,surv)
  rmst <- sum(area)
  
  n.risk <- ft$n.risk[ind]
  n.event <- ft$n.event[ind]
  
  wk.var <- c(if_else((n.risk-n.event)==0, 0, n.event /(n.risk *(n.risk - n.event))),0)
  rmst.var <- sum(cumsum(rev(area[-1]))^2 * rev(wk.var)[-1])
  
  return(list(rmst=rmst, rmst.var=rmst.var))
}

RMST_comparison <- function(data, tau, k_tau, family=NULL, method='indep', alpha=0.05, ensemble=F, weight=c(1,1,1), tol=1e-6){
  
  dat1 <- data %>% filter(trt==1)
  dat0 <- data %>% filter(trt==0)
  
  if(method=='indep'){
    rmst_res1 <- RMST_indep(data=dat1, tau=tau)
    rmst1 <- rmst_res1$rmst
    rmst.var1 <- rmst_res1$rmst.var
    rmst_res0 <- RMST_indep(data=dat0, tau=tau)
    rmst0 <- rmst_res0$rmst
    rmst.var0 <- rmst_res0$rmst.var
  }else if(method!='indep'&ensemble==F){
    rmst1 <- RMST_cpp(data=dat1, tau=tau, method=method, theta=BiCopTau2Par(family=family, tau=k_tau), family=family)
    rmst.var1 <- RMST_var_cpp(data=dat1, tau=tau, method=method, theta=BiCopTau2Par(family=family, tau=k_tau), family=family)
    rmst0 <- RMST_cpp(data=dat0, tau=tau, method=method, theta=BiCopTau2Par(family=family, tau=k_tau), family=family)
    rmst.var0 <- RMST_var_cpp(data=dat0, tau=tau, method=method, theta=BiCopTau2Par(family=family, tau=k_tau), family=family)
  }else if(method!='indep'&ensemble==T){
    rmst1 <- RMST_ENS_cpp(data=dat1, tau=tau, method=method, weight=weight)
    rmst.var1 <- RMST_var_cpp(data=dat1, tau=tau, method=method, ensemble=T, weight=weight)
    rmst0 <- RMST_ENS_cpp(data=dat0, tau=tau, method=method, weight=weight)
    rmst.var0 <- RMST_var_cpp(data=dat0, tau=tau, method=method, ensemble=T, weight=weight)
  }else{
    stop('Check method')
  }

  rmst.diff <- rmst1 - rmst0
  rmst.diff.se <- sqrt(rmst.var1 + rmst.var0)
  rmst.diff.low <- rmst.diff - qnorm(1 - alpha/2) * rmst.diff.se
  rmst.diff.upp <- rmst.diff + qnorm(1 - alpha/2) * rmst.diff.se
  rmst.diff.pval <- pnorm(-abs(rmst.diff)/rmst.diff.se) * 2
  rmst.diff.result <- c(rmst.diff, rmst.diff.low, rmst.diff.upp, rmst.diff.pval)
  out <- rbind(rmst.diff.result)
  rownames(out) <- c("RMST (arm=1)-(arm=0)")
  colnames(out) <- c("Est.", paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""), 
                     paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""), "p")
  out
}



