# MICROBENCHMARK TESTS between bonorico/gcipdr and bonofi/gcipdr_test@optimization_Rccp

# install dependencies

url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
pkgFile <- "JohnsonDistribution_0.24.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)

library(microbenchmark)
library(remotes)

pak::pak("bonorico/gcipdr")

##### chnaging session might disrupt path to Rtools
Sys.setenv(PATH = paste(
  "C:/RBuildTools/4.4/usr/bin",
  "C:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/bin",
  Sys.getenv("PATH"),
  sep = ";"
))

# Verify make is found
Sys.which("make")
Sys.which("gcc")

Rcpp::compileAttributes()
devtools::load_all()

#pak::pak("bonofi/gcipdr_test@optimization_Rccp2")

library(gcipdr)
library(gcipdrtest)



#### START TESTS

testdat <- mtcars[, 1:4]

# Compare if Rccp fine tuning 
# speeds up routine  
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
