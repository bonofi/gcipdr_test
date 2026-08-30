# MICROBENCHMARK TESTS between bonorico/gcipdr and bonofi/gcipdr_test@optimization_kruskal_inits

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
pak::pak("bonofi/gcipdr_test@optimization_kruskal_inits")

library(gcipdr)
library(gcipdrtest)



#### START TESTS

testdat <- mtcars[, 1:2]
apply(testdat, 2, mean)
apply(testdat, 2, sd)

test <- gcipdrtest::Simulate.data.given.IPD(
  testdat, H=5, stochastic.integration = TRUE, 
  SI_k = 50000, method = 10, # method 10 user-defined marginals 
  checkdata = TRUE,
  tabulate.similar.data = TRUE,
  user_defined_marginals = list(
    function(x) qnorm(x, mean = 20, sd = 6),
    function(x) qlnorm(x, meanlog = log(6), sdlog = log(1.8))
  )
  
)


Jparms1 <- JohnsonFit(testdat$mpg)



test2 <- gcipdrtest::Simulate.data.given.IPD(
  testdat, H=5, stochastic.integration = TRUE, 
  SI_k = 50000, method = 10, # method 10 user-defined marginals 
  checkdata = TRUE,
  tabulate.similar.data = TRUE,
  user_defined_marginals = list(
    function(x) qJohnson(x, parms = Jparms1),
    function(x) qnorm(x, mean = 20, sd = 6)
    )
  )






res <- microbenchmark::microbenchmark(
  {
    set.seed(608, "L'Ecuyer")
    gcipdr::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, 
                                    SI_k = 50000, method = 3, checkdata = TRUE,
                                    tabulate.similar.data = TRUE)},
  {
    set.seed(608, "L'Ecuyer")
    gcipdrtest::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, 
                                        SI_k = 50000, method = 3, checkdata = TRUE,
                                        tabulate.similar.data = TRUE)},
  times = 30L,
  check = "equal"
)

print(res)

boxplot(res, names = c("gcipdr", "gcipdrtest"))