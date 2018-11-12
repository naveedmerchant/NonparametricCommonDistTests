#install.packages("rmutil")
library(rmutil)
source("MarginalLikIntfunctions.R")
set.seed(10000)
dlength <- 1000
LBF <- c()
LBF2 <- c()
dataset1 <- rcauchy(dlength)
dataset2 <- rcauchy(dlength)
for(i in 1:30)
{
  dataset1 <- sample(dataset1)
  dataset2 <- sample(dataset2)
  XT1 <- dataset1[1:(dlength*.3)]
  XV1 <- dataset1[-(1:dlength*.3)]
  XT2 <- dataset2[1:(dlength*.3)]
  XV2 <- dataset2[-(1:dlength*.3)]
  ExpectedKernML <- logmarg.kernMCimport(XT1,XV1,iter = 1, importsize = 100) + logmarg.kernMCimport(XT2,XV2,iter = 1,importsize = 100)
  ExpectedKernML2 <- logmarg.kernMCimport(c(XT1,XT2),c(XV1,XV2),iter = 1,importsize = 100)
  LBF[i] <- ExpectedKernML - ExpectedKernML2
}

plot(LBF)

#Pareto distribution has finite mean but infinite variance, if dispersion is bigger than 1 but less than 2.
#Thus...

var(rpareto(1000,m = 1.5, s = 1.5))



set.seed(10000)
dlength <- 100
LBF <- c()
for(i in 1:100)
{
  dataset1 <- rpareto(dlength, m=1.5, s=1.5)
  dataset2 <- rpareto(dlength, m=1.5, s=1.5)
  
  XT1 <- dataset1[1:(dlength*.3)]
  XV1 <- dataset1[-(1:dlength*.3)]
  XT2 <- dataset2[1:(dlength*.3)]
  XV2 <- dataset2[-(1:dlength*.3)]
  ExpectedKernML <- logmarg.kern(XT1,XV1)[[2]] + logmarg.kern(XT2,XV2)[[2]]
  ExpectedKernML2 <- logmarg.kern(c(XT1,XT2),c(XV1,XV2))[[2]]
  LBF[i] <- ExpectedKernML - ExpectedKernML2
}

plot(LBF)

#Forcing finite variance and mean but infinite third moment

set.seed(10000)
dlength <- 100
LBF <- c()
for(i in 1:100)
{
  dataset1 <- rpareto(dlength, m=2.5, s=2.5)
  dataset2 <- rpareto(dlength, m=2.5, s=2.5)
  
  XT1 <- dataset1[1:(dlength*.3)]
  XV1 <- dataset1[-(1:dlength*.3)]
  XT2 <- dataset2[1:(dlength*.3)]
  XV2 <- dataset2[-(1:dlength*.3)]
  ExpectedKernML <- logmarg.kern(XT1,XV1)[[2]] + logmarg.kern(XT2,XV2)[[2]]
  ExpectedKernML2 <- logmarg.kern(c(XT1,XT2),c(XV1,XV2))[[2]]
  LBF[i] <- ExpectedKernML - ExpectedKernML2
}

plot(LBF)

#similar problems persist


