#' This is a function provided for simulation
#' 
#' @param N sample size. default value is 5000.
#' @import mvtnorm

simul<-function(N = 5000){
  #generate z
  p=5
  Sigma_z1<-AR1(0.6,p)
  
  z= rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1)
  
  
  ##generate delta and time
  U=runif(N, 0,1)
  
  pre_time=rep(0, N)
  for (i in 1:(N)) {
    f=function(t) {
      integrand <- function(x) {0.5*exp(gamma_subject[i]+z[i,1]-z[i,3]+sin(3*pi*x/4)*(x<3)*z[i,2]+z[i,4]-z[i,5])}
      Lambda=integrate(integrand, lower = 0, upper = t)$value
      Lambda+log(1-U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f,  lower = 0, upper = 4), silent=TRUE))
    if (class(r1) == "try-error"){    
      pre_time[i]=4
    }
    else pre_time[i]=uniroot(f,  lower = 0, upper = 4)$root
  }
  
  pre_censoring=runif(N,0,3)
  pre_censoring=pre_censoring*(pre_censoring<3)+3*(pre_censoring>=3)
  tcens=(pre_censoring<pre_time) # censoring indicator
  delta=1-tcens
  time=pre_time*(delta==1)+pre_censoring*(delta==0)
  
  return(list(z=z,delta=delta,time=time))
  
  
}