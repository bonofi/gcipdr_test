# MICROBENCHMARK TESTS between bonorico/gcipdr and bonofi/gcipdr_test

# install dependencies

url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
pkgFile <- "JohnsonDistribution_0.24.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)

library(tidyverse)
library(microbenchmark)
library(remotes)
library(SuppDists)
SuppDists::moments(testdat$mpg)

pak::pak("bonorico/gcipdr")

library(gcipdr)
library(gcipdrtest)



#### START TESTS

testdat <- rbind(
  mtcars[, 1:4],
  mtcars[, 1:4],
  mtcars[, 1:4],
  mtcars[, 1:4]
)

res <- microbenchmark::microbenchmark(
  {
    gcipdr::Simulate.data.given.IPD(
      testdat, H=5000, method = 3)},
  {
    gcipdrtest::Simulate.data.given.IPD(
      testdat, H=5000, method = 3
                                        )},
  times = 30L
)


print(res)
boxplot(res, names = c("gcipdr", "gcipdrtest"))