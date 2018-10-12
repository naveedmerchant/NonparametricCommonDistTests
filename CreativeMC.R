dataset<- rnorm(500)
X1 <- dataset[1:(length(dataset)*.3)]
X2 <- dataset[-(1:(length(dataset)*.3))]
B <- unname((quantile(X2,probs = .75) - quantile(X2,probs = .25))/1.35)
k <- length(X1)
n <- length(X2)
#Normconst <- 10
sum2 <- 0
iter <- 10000
for(i in 1:iter)
{
  sum1 <- 0
  for(j in 1:n)
  {
    l <- sample(1:k,1)
    sum1 <- sum1 +  (X2[j] - X1[l])^2
  }
  sum1 <- lgamma((n+3)/2) - ((n+3)/2)*log(.5*(2*B^2 + sum1)) 
  sum2 <- sum2 + exp(sum1)
}
logMCinteg <-log(B / sqrt(pi)) -(n/2)*log(2*pi) + log(sum2) -log(iter)
#Split factorial into components
#Split exponentiation into components


#Compare with importance sampling estimate:


#Computes prior after adjusting for drawing from cauchy

#The next few lines evaluate the likelihood of the model under the gaussian kernel.

Loglist <- c()

for(G in 1:100)
{
  cauchsamp<- rcauchy(1000)

  poscauchsamp <- cauchsamp[cauchsamp > 0] 

  importancepart <- ((2*B*(1/poscauchsamp^2)*(exp(-(B^2)/(poscauchsamp)^2)) / sqrt(pi))) / (2*dcauchy(poscauchsamp))


  sum <- 0
  prod <- 0
  sum2 <- 0
  prod2 <- 0

  prodlist1<-c()
  prodlist2<-c()

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
      prod <- prod + log(sum) + log((1/sqrt(2*pi))) + log((k*poscauchsamp[z])^(-1))
    }
    prodlist1[z] <- prod + log(importancepart[z])
  }
  log(sum(exp(prodlist1)))
}