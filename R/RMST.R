#' @export
#'
#' @import survival
#' @import dplyr
#' @import VineCopula
#' @import copula
#' @importFrom foreach %dopar% foreach
#' @importFrom iterators icount
#' @importFrom doParallel registerDoParallel

RMST_var <- function(data, method, tau, theta=1, family=3, ensemble=FALSE, theta_vec=c(2/3,2,6), weight=c(1,1,1), tol=1e-6,
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
                                .combine = "c",
                                .packages = "RMSTdepC") %dopar% {
    boot_data = data[sample(1:nrow(data), replace=T),]
    if(ensemble==FALSE){
      result = RMST_cpp(boot_data, method, tau, family, theta, tol);
    }else{
      result = RMST_ENS_cpp(boot_data, method, tau, theta_vec, weight, tol);
    }
    result;
  }
  var(output)
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

RMST_comparison <- function(data, tau, k_tau, family=3, method='indep', alpha=0.05, ensemble=F, theta_vec=c(2/3,2,6), weight=c(1,1,1), tol=1e-6, parallel=F,
                            n_boots = 999, n_cores = parallel::detectCores()){
  
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
    rmst1 <- RMST_cpp(data=dat1, tau=tau, method=method, family=family, theta=BiCopTau2Par(family=family, tau=k_tau), tol=tol)
    rmst0 <- RMST_cpp(data=dat0, tau=tau, method=method, family=family, theta=BiCopTau2Par(family=family, tau=k_tau), tol=tol)

    if(parallel){
      rmst.var1 <- RMST_var(data=dat1, tau=tau, method=method, family=family, theta=BiCopTau2Par(family=family, tau=k_tau), tol=tol, n_boots=n_boots, n_cores=n_cores)
      rmst.var0 <- RMST_var(data=dat0, tau=tau, method=method, family=family, theta=BiCopTau2Par(family=family, tau=k_tau), tol=tol, n_boots=n_boots, n_cores=n_cores)
    }else{
      rmst.var1 <- RMST_var_cpp(data=dat1, tau=tau, method=method, family=family, theta=BiCopTau2Par(family=family, tau=k_tau), tol=tol, n_boots=n_boots)
      rmst.var0 <- RMST_var_cpp(data=dat0, tau=tau, method=method, family=family, theta=BiCopTau2Par(family=family, tau=k_tau), tol=tol, n_boots=n_boots)
    }
      
  }else if(method!='indep'&ensemble==T){
    rmst1 <- RMST_ENS_cpp(data=dat1, tau=tau, method=method, theta_vec=theta_vec, weight=weight, tol=tol)
    rmst0 <- RMST_ENS_cpp(data=dat0, tau=tau, method=method, theta_vec=theta_vec, weight=weight, tol=tol)
    
    if(parallel){
      rmst.var1 <- RMST_var(data=dat1, tau=tau, method=method, ensemble=T, theta_vec=theta_vec, weight=weight, tol=tol, n_boots=n_boots, n_cores=n_cores)
      rmst.var0 <- RMST_var(data=dat0, tau=tau, method=method, ensemble=T, theta_vec=theta_vec, weight=weight, tol=tol, n_boots=n_boots, n_cores=n_cores)
    }else{
      rmst.var1 <- RMST_var_cpp(data=dat1, tau=tau, method=method, ensemble=T, theta_vec=theta_vec, weight=weight, tol=tol, n_boots=n_boots)
      rmst.var0 <- RMST_var_cpp(data=dat0, tau=tau, method=method, ensemble=T, theta_vec=theta_vec, weight=weight, tol=tol, n_boots=n_boots)
    }

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



