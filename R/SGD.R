#' This is the basic SGD version.
#' 
#' @param delta event indicator
#' @param z Covariate matrix
#' @param time observed event time
#' @param knot number of basis functions for time-varing effects. the default value is 10.
#' @param facility strata. 
#' @param M_stop maximum stopping iterations
#' @param tol the convergence threshold
#' 
#' @import splines
#' @export 

SGD = function(delta, z, time, knot = 10 ,facility = NULL, theta_SGD = NULL, M_stop = 10000, tol= 10^(-6)  ) {
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

  
  data=cbind(time,delta,facility, bs8,z)
  data=data[order(facility),]
  m=as.vector(sapply(split(data[,3],factor(data[,3])),length))
  
  m_stratify=0
  theta_stratify=matrix( rep(0, knot*p), nrow=p)  #dim P*knot
  theta_SGD_all=NULL
  likelihood_SGD_all=dloglik_likelihood_gradient(knot,facility,delta,z,bs8,theta_stratify,N)
  
  
  stage2_key=FALSE
  while (stage2_key==FALSE){
    m_stratify=m_stratify+1
    
    set.seed(m_stratify)
    bootsample3<-stra_sampling(m,10)
    data_sub4=data[bootsample3,]
    
    time_sub4=data_sub4[,1]
    delta_sub4=data_sub4[,2]
    facility_sub4=data_sub4[,3]
    b_spline_sub4=data_sub4[,4:(knot+3)]
    z_sub4=data_sub4[,-(1:(knot+3))]
    
    
    #note that here do the sort again.
    delta_sub4 = delta_sub4[order(time_sub4)]
    facility_sub4=facility_sub4[order(time_sub4)]
    z_sub4 = z_sub4[order(time_sub4),]
    b_spline_sub4=b_spline_sub4[order(time_sub4),]
    time_sub4 = time_sub4[order(time_sub4)]
    n_sub4=length(facility_sub4)
    
    result=GDboost_gradient(knot,rate,facility_sub4,delta_sub4,z_sub4,b_spline_sub4,theta_stratify)
    
    theta_stratify=theta_stratify-rate*matrix(result$L1,nrow=p,byrow = TRUE)
    
    theta_SGD_all[[m_stratify]]=theta_stratify
    
    likelihood_stratify=dloglik_likelihood_stratify(knot,facility,delta,z,bs8,theta_stratify,N)
    likelihood_SGD_all=c(likelihood_SGD_all, likelihood_stratify)
    
    
    beta_SGD=matrix( rep(0, N*p), nrow=N)
    
    
    for (j in 1:p) {
      beta_SGD[,j]=bs8%*%theta_stratify[j,]
      
    }
    
    
    if (m_stratify>=(10+track)) {
      llk.diff.all=NULL
      for (key in 1:track) {
        llk.diff.all = c(llk.diff.all, likelihood_SGD_all[m_stratify-key+1]-likelihood_SGD_all[m_stratify-key])
      }
      llk.diff=max(llk.diff.all)
      llk.diff_null = likelihood_SGD_all[m_stratify]-likelihood_SGD_all[1]
      dist=theta_SGD_all[[m_stratify]]-theta_SGD_all[[m_stratify-1]]
      if(abs(llk.diff/llk.diff_null) < tol) {
        stage2_key=TRUE 
        break
      }
    }
    
    if (m_stratify==M_stop){
      stage2_key=TRUE
      break
    }
  } #end while
  
  theta.SGD.all=rbind(theta.SGD.all, theta_stratify)
  
  temp=ddloglik(knot,facility,delta,z,bs8,theta_stratify,length(unique(facility)))
  
  test_SGD=rep(0,p)
  constrain=-diff(diag(knot*1),differences=1)
  j=0
  repeat{
    j=j+1
    
    theta_SGD_j= theta_stratify[j,]
    
    L2=solve(temp$GVG)[((j-1)*knot+1):((j-1)*knot+knot),((j-1)*knot+1):((j-1)*knot+knot)]
    test_contrast=t(constrain%*%theta_SGD_j)%*%solve(constrain%*%L2%*%t(constrain))%*%(constrain%*%theta_SGD_j)
    test_SGD[j]=1-pchisq(test_contrast, (knot-1))
    if(j==p) break
  }
  
  test_SGD_all=rbind(test_SGD_all, test_SGD)
  
  result = list(theta = theta.SGD.all, test = test_SGD_all)
  
  return(result)
  
}

