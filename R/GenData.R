#' @export
#'
#' @import copula
#' @import VineCopula
#' @import dplyr

time_gen <- function(n, dist, param){
  if(dist=='exp'){
    rexp(n, rate=param)
  }else if(dist=='weibull'){
    rweibull(n, shape=param[1], scale=param[2])
  }
}

indep_gen_data <- function(n, dist_1=c('exp','exp'), dist_2=c('exp','exp'),
                           param_1=list(0.2, 0.2),
                           param_2=list(0.2, 0.2)){
  
  t1 <- time_gen(n, dist_1[1], param_1[[1]]); c1 <- time_gen(n, dist_1[2], param_1[[2]]) # rweibull(num_obs, shape = 3, scale = 9)
  t2 <- time_gen(n, dist_2[1], param_2[[1]]); c2 <- time_gen(n, dist_2[2], param_2[[2]]) # rweibull(num_obs, shape = 0.5, scale = 20)
  data <- data.frame(t = c(t1, t2), 
                     c = c(c1, c2), 
                     trt = c(rep(1, n), rep(0, n))) %>% 
    mutate(time = pmin(t, c), status = ifelse(t <= c, 1, 0)) %>%
    select(time, status, trt)
  
  return(data)
}

dep_gen_data <- function(n, dist_1=c('exp','exp'), dist_2=c('exp','exp'),
                         param_1=list(list(rate=0.2), list(rate=0.2)),
                         param_2=list(list(rate=0.2), list(rate=0.2)), family, tau_given){
  
  copula_type <- switch(family-2, 'clayton', 'gumbel', 'frank')
  
  theta_given <- BiCopTau2Par(family = family, tau = tau_given)
  cop <- archmCopula(copula_type, param = theta_given)
  mv1 <- mvdc(cop, margins = dist_1, paramMargins = param_1)
  mv2 <- mvdc(cop, margins = dist_2, paramMargins = param_2)
  dep_tc1 <- rMvdc(n, mv1); dep_tc2 <- rMvdc(n, mv2)
  data <- data.frame(t = c(dep_tc1[,1], dep_tc2[,1]), 
                     c = c(dep_tc1[,2], dep_tc2[,2]), 
                     trt = c(rep(1, n), rep(0, n))) %>% 
    mutate(time = pmin(t, c), status = ifelse(t <= c, 1, 0)) %>%
    select(time, status, trt, t, c)
  
  return(data)
}

