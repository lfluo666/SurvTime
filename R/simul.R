simul <- function(N = NULL){
  #generate z
  if(is.null(N)){
    N=5000
  }
  p=5
  Sigma_z1<-AR1(0.6,p)
  
  z= rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1)
  
  F=1 ###number of facility
  n_f = rep(N, F)#rpois(F, lambda = 1000)  #sample size for each facility
  gamma = rnorm(F, mean=0, sd=0.5)
  range(gamma)
  gamma_subject=rep(gamma,n_f)
  
  ##generate delta and time
  U=runif(N, 0,1)
  
  pre_time=rep(0, N)
  for (i in 1:(N)) {
    f=function(t) {
      #integrand <- function(x) {0.5*exp(gamma_subject[i]+z[i,1]-z[i,3]+sin(3*pi*x/4)*(x<3)*z[i,2]+exp(-x)*(x<3)*z[i,4]+z[i,5])}
      integrand <- function(x) {0.5*exp(gamma_subject[i]+z[i,1]-z[i,3]+sin(3*pi*x/4)*(x<3)*z[i,2]-z[i,4]+z[i,5])}
      
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

#' AR1
AR1 <- function(tau, m) {
  if(m==1) {R <- 1}
  if(m > 1) {
    R <- diag(1, m)
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        R[i,j] <- R[j,i] <- tau^(abs(i-j))
      }
    }
  }
  return(R)
}

#' stra_sampling
stra_sampling<-function(size, prop){
  m<-length(size)
  mark<-cumsum(size)
  total<-sum(size)
  sample.id<-NULL
  sample.id<-c(sample.id,sample.int(size[1],size=floor(size[1]/prop),replace=FALSE))
  
  if (m>1) {
    for(i in 2:m){
      sample.id<-c(sample.id,mark[i-1]+sample.int(size[i],size=floor(size[i]/prop),replace=FALSE))
    }
  }
  return(sample.id)
}


