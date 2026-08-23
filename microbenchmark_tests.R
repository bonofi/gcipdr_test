# MICROBENCHMARK TESTS between bonorico/gcipdr and bonofi/gcipdr_test@optimization_test

# install dependencies

url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
pkgFile <- "JohnsonDistribution_0.24.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)

library(microbenchmark)

pak::pak("bonorico/gcipdr")
install_github("bonofi/gcipdr_test", ref="optimization_test")

library(gcipdr)
library(gcipdrtest)



#### START TESTS

testdat <- mtcars[, 1:4]

# Compare Base R vs dplyr methods
res <- microbenchmark(
  gcipdr::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, SI_k = 50000, method = 3),
  gcipdrtest::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, SI_k = 50000, method = 3),
  times = 10L,
  check = "equal"
)


  
boxplot(res, names = c("gcipdr", "gcipdrtest"))


