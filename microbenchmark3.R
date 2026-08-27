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
install_github("bonofi/gcipdr_test", ref="optimization_Rccp")

library(gcipdr)
library(gcipdrtest)



#### START TESTS

testdat <- mtcars[, 1:4]

# Compare if unification of rmvnorm generation in stochastic integration 
# speeds up routine (expected x2) 
# (commit: 1c26551488b77a6b5f6e6facc89bd07fb59253df)

custom_check <- function(res){
  
  tol <- 0.2 # tolerance parameter
  x <- res[[1]]
  y <- res[[2]]
  
  gc_param_x <- res[[1]]$copula.parameters
  gc_param_y <- res[[2]]$copula.parameters
  
  all(gc_param_x[lower.tri(gc_param_x)] - gc_param_y[lower.tri(gc_param_y)])  
}


res <- microbenchmark(
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
  #check = custom_check
)

print(res)

boxplot(res, names = c("gcipdr", "gcipdrtest"))
