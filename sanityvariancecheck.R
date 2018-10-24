source("MarginalLikIntfunctions.R")
set.seed(10000)
dlength <- 100
LBF <- c()
for(i in 1:100)
{
  dataset1 <- rcauchy(dlength)
  dataset2 <- rcauchy(dlength)
  
  XT1 <- dataset1[1:(dlength*.3)]
  XV1 <- dataset1[-(1:dlength*.3)]
  XT2 <- dataset2[1:(dlength*.3)]
  XV2 <- dataset2[-(1:dlength*.3)]
  ExpectedKernML <- logmarg.kern(XT1,XV1)[[2]] + logmarg.kern(XT2,XV2)[[2]]
  ExpectedKernML2 <- logmarg.kern(c(XT1,XT2),c(XV1,XV2))[[2]]
  LBF[i] <- ExpectedKernML - ExpectedKernML2
}

plot(LBF)

