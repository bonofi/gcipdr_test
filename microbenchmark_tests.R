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

# Compare if vectorization in stochastic integration has speed benefits 
# (commit: 32e7fae2e3998d5366de646c8468ce9c244853b5)

res <- microbenchmark(
  {
    set.seed(608, "L'Ecuyer")
    gcipdr::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, SI_k = 50000, method = 3)},
  {
    set.seed(608, "L'Ecuyer")
    gcipdrtest::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, SI_k = 50000, method = 3)},
  times = 30L,
  check = "equal"
)

print(res)
  
boxplot(res, names = c("gcipdr", "gcipdrtest"))


