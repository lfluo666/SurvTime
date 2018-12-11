NR = function(delta, z, knot = 10, facility = NULL, theta_NR = NULL, M_stop = 1000) {
  
  K=knot
  N=nrow(z)
  p=ncol(z)
  theta.SGD.all=NULL
  test_SGD_all=NULL
  rate=0.001 
  track=5
  if(is.null(facility)){
    F=1 ##set facility = 1
    n_f = rep(N,F) ## sample size for each facility
    F_pre=1:F 
    facility=rep(F_pre, n_f)
  }
  ##(1) when they give facility, how to determine F,n_f,???
  
  
  #compute bs8
  delta = delta[order(time)]
  facility=facility[order(time)]
  z = z[order(time),]
  time = time[order(time)]
  
  
  time2=time[delta==1]
  knot_set=quantile(time2,prob=seq(1:(knot-4))/(knot-3))
  bs7=splines::bs(time,df=knot, knot=knot_set, intercept=TRUE, degree=3)
  bs8=matrix(bs7, nrow=N) # for build beta_t
  
  b_spline=bs8

  if(is.null(theta_NR)){
    theta_NR = matrix(0,nrow = p,ncol = knot)
  }
  
  number_facility = length(unique(facility))
  result1=NRrepeat(N,knot , facility,delta, z, bs8, theta_NR, number_facility ,M_stop)
  
  # two results to return (3 is temp.GVG to compute test_NR_all)
  theta_NR = result1["theta"]
  likelihood_NR_all = result1["likelihood"]
  
  theta.NR.all=rbind(theta.NR.all, theta_NR)
  
  return(theta_NR)
}