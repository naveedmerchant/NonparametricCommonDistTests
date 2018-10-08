#Run this as a loop and examine shape of density of BF's
LBF<-c()
set.seed(900)
for(D in 1:100)
{
  dataset<- rnorm(500)
  X1 <- dataset[1:(length(dataset)*.3)]
  X2 <- dataset[-(1:(length(dataset)*.3))]
  datasetperm <- sample(dataset)
  X1perm <- datasetperm[1:(length(dataset)*.3)]
  X2perm <- datasetperm[-(1:(length(dataset)*.3))]

  cauchsamp<- rcauchy(1000)

  poscauchsamp <- cauchsamp[cauchsamp > 0] 

#We draw from truncated cauchy, this is because we want our bandwidth to be always positive. 

B <- (quantile(X2,probs = .75) - quantile(X2,probs = .25))/1.35

#Used in prior for bandwidth

importancepart <- ((2*B*(1/poscauchsamp^2)*(exp(-(B^2)/(poscauchsamp)^2)) / sqrt(pi))) / (2*dcauchy(poscauchsamp))

#Computes prior after adjusting for drawing from cauchy

#The next few lines evaluate the likelihood of the model under the gaussian kernel.


k <- length(X1)
n <- length(X2)
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
for(z in 1:length(poscauchsamp))
{
  prod <- 0
  for(j in 1:n)
  {
    sum <- 0
    for(i in 1:k)
    {
      sum <- sum + exp(-.5*((X2perm[j]-X1perm[i])/poscauchsamp[z])^2)
    }
    prod <- prod + log(sum) + log((1/sqrt(2*pi))) + log((k*poscauchsamp[z])^(-1))
  }
  prodlist2[z] <- prod + log(importancepart[z])
}
BF1 <- sum(exp(prodlist1-max(prodlist2))) / sum(exp(prodlist2-max(prodlist2)))
logBF1 <- log(BF1)
LBF[D]<-logBF1
}

summary(LBF)
sd(LBF)
median(LBF)
mean(LBF)
#This is a little more numerically stable.

#Epanech still does far worse than Gaussian Kernel

#The bayes factor isn't 0, but the marginal likelihood of the epanech kernel is 
#Too small too compute, the computer simply guesses it to be 0 unless it is multiplied by a big number due to underflow

#If n is low then it seems to be that the gaussian kernel is favored.
#If n is high it seems to be that the epanech kernel is favored

#This calculation is very time expensive!

#One of the permutations is favored quite harshly!

#We repeat this computation many times and observe the results
