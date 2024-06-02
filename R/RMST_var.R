#' @export
#'
#' @importFrom foreach %dopar% foreach
#' @importFrom iterators icount
#' @importFrom doParallel registerDoParallel
#' @importFrom stats rnorm

RMST_var <- function(data, method, tau, theta, family,
                     n_boots = 1000,
                     n_cores = parallel::detectCores()) {
  
  # Construct cluster
  cl = parallel::makeCluster(n_cores)
  
  # After the function is run, close the cluster.
  on.exit(parallel::stopCluster(cl))
  
  # Register parallel backend
  doParallel::registerDoParallel(cl)
  
  
  # Compute estimates
  estimates <- foreach::foreach(i = iterators::icount(n_boots), # Perform n simulations
                                .combine = "rbind",           # Combine results
                                # Self-load
                                .packages = "Rcpp2doParallel") %dopar% {
                                  boot_data = data[sample(1:nrow(data), replace=T),]
                                  result = RMST_cpp(boot_data, method, tau, theta, family);
                                  result
                                }
  
  estimates
}

RMST <- function(data, method, tau, theta, family) {
  RMST_cpp(boot_data, method, tau, theta, family)
}

