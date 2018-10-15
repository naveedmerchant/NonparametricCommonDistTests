logmarg.kern=function(y,x,prior=1){
  
  out=optim(.4,loglike.KGauss,method="L-BFGS-B",lower=.0001,upper=5,y=y,x=x)
  
  h=out$par
  
  cons=-out$val
  
  stat=integrate(integrand.Gauss,lower=0.0001,upper=5,y=y,x=x,cons=cons,prior=prior)$val
  
  list(h,cons+log(stat))
  
}


loglike.KGauss=function(h,y,x){
  
  n=length(x)
  
  m=length(y)
  
  nh=length(h)
  
  M=t(matrix(y,m,n))
  
  llike=1:nh
  
  for(j in 1:nh){
    
    M1=(x-M)/h[j]
    
    M1=dnorm(M1)/h[j]
    
    fhat=as.vector(M1 %*% matrix(1,m,1))/m
    
    fhat[fhat<10^(-320)]=10^(-320)
    
    llike[j]=sum(log(fhat))
    
  }
  
  -llike
  
}



integrand.Gauss=function(h,y,x,cons,prior){
  
  n=length(x)
  
  R=quantile(y,probs=c(.25,.75))
  
  R=R[2]-R[1]
  
  beta=R/1.35
  
  beta1=beta*log(2)/sqrt(qgamma(.5,.5,1))
  
  Prior=beta1*exp(-beta1/h)/h^2
  
  if(prior==1) Prior=(2*beta/sqrt(pi))*h^(-2)*exp(-beta^2/h^2)
  
  arg=-loglike.KGauss(h,y,x)-cons
  
  arg[arg>700]=700
  
  f=exp(arg)*Prior
  
  f
  
}