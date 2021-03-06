library(matrixStats)
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

logmarg.kernMC=function(X1,X2,iter = 10000)
{
  require(matrixStats)
  n <- length(X2)
  k <- length(X1)
  sum1 <- c()  
  sum1list <- c()
  R = quantile(X1,probs=c(.25,.75))
  R=R[2]-R[1]
  R = unname(R)
  B = R / 1.35
  #browser()
  for(i in 1:iter)
  {
    l <- sample(1:k,n, replace = TRUE)
    sum1 <- sum((X2 - X1[l])^2)
    sum1 <- lgamma((n-1)/2) - ((n-1)/2)*log(.5*(2*B^2 + sum1)) 
    sum1list[i] <- sum1
  }
  logMCinteg <-log(B / sqrt(pi)) - (n/2)*log(2*pi) + logSumExp(sum1list) - log(iter)
  return(logMCinteg)
}


logmarg.kernMCimport=function(X1,X2,iter = 50,importsize = 200)
{
  R = quantile(X1,probs=c(.25,.75))
  R=R[2]-R[1]
  R = unname(R)
  B = R / 1.35
  Loglist <- c()
  k <- length(X1)
  n <- length(X2)
  for(G in 1:iter)
  {
    
    cauchsamp<- rcauchy(importsize)
    
    poscauchsamp <- cauchsamp[cauchsamp > 0] 
    
    importancepart <- ((2*B*(1/poscauchsamp^2)*(exp(-(B^2)/(poscauchsamp)^2)) / sqrt(pi))) / (2*dcauchy(poscauchsamp))
    
    prodlist1<-c()
    
    for(z in 1:length(poscauchsamp))
    {
      prod <- 0
      for(j in 1:n)
      {
        sum <- 0
        for(i in 1:k)
        {
          sum <- sum + exp(-.5*((X2[j]-X1[i])/poscauchsamp[z])^2)
        }
        prod <- prod + log(sum) + log((1/sqrt(2*pi))) - log((k*poscauchsamp[z]))
      }
      prodlist1[z] <- prod + log(importancepart[z])
    }
    Loglist[G] <- logSumExp(prodlist1)
  }
  return(Loglist)
}
