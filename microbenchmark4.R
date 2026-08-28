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
# (commit: e2261b14c89160b133649fd59b77ddf85cd0b6cb)

custom_check <- function(res){
  
  tol <- 0.01 # tolerance parameter for correlations
  x <- res[[1]]
  y <- res[[2]]
  
  browser()
  
  gc_param_x <- res[[1]]$copula.parameters
  gc_param_y <- res[[2]]$copula.parameters
  
  corr_x <- res[[1]]$is.data.similar$lower.triangular.Rx$diff
  corr_y <- res[[2]]$is.data.similar$lower.triangular.Rx$diff
  
  # check copula parameters difference with tolerance value
  check1 <- all(
    (gc_param_x[lower.tri(gc_param_x)] - gc_param_y[lower.tri(gc_param_y)]) <= tol
    )  
  
  # # check correlation differences with tolerance value
  # check2 <- all(
  #   (corr_x - corr_y) <= tol
  # )  
  # 
  # check1 & check2

  check1
  }
#
#
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
  times = 30L
  #check = "equal"
)

print(res)

boxplot(res, names = c("gcipdr", "gcipdrtest"))


### check equivalence with custom check (increase H) on 1-time benchmark
for ( i in 1:10){
  
  set.seed(608+i, "L'Ecuyer")
  ref <- gcipdr::Simulate.data.given.IPD(testdat, H=1, stochastic.integration = TRUE, 
                                         SI_k = 50000, method = 3, checkdata = TRUE,
                                         tabulate.similar.data = TRUE)
  
  set.seed(608+i, "L'Ecuyer")
  new <- gcipdrtest::Simulate.data.given.IPD(testdat, H=1, stochastic.integration = TRUE, 
                                             SI_k = 50000, method = 3, checkdata = TRUE,
                                             tabulate.similar.data = TRUE)
  
  custom_check(list(ref, new)) 
}