# MICROBENCHMARK TESTS between bonorico/gcipdr and bonofi/gcipdr_test@optimization_kruskal_inits

# install dependencies

url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
pkgFile <- "JohnsonDistribution_0.24.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)

library(microbenchmark)
library(remotes)

pak::pak("bonorico/gcipdr")
pak::pak("bonofi/gcipdr_test@optimization_kruskal_inits")

library(gcipdr)
library(gcipdrtest)



#### START TESTS

testdat <- mtcars[, 1:5]

# Compare if Kruskal initialization  
# speeds up routine (expected x2) 
# (commit: c69fa3257b6076db80cfc33868b14938aa395715)




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