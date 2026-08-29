# MICROBENCHMARK TESTS between bonorico/gcipdr and bonofi/gcipdr_test@optimization_parallel

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

install_github("bonofi/gcipdr_test", ref="optimization_parallel")

library(gcipdr)
library(gcipdrtest)
library(future)
library(future.apply)

#### START TESTS

testdat <- mtcars[, 1:4]

# Compare if future parallelization 
# speeds up routine (expected x2) 
# (commit: 407b55895dba452f6de14c0f641eb36233db9cc7)
# 
res <- microbenchmark::microbenchmark(
  {
    set.seed(608, "L'Ecuyer")
    gcipdr::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, 
                                    SI_k = 50000, method = 3, checkdata = TRUE,
                                    tabulate.similar.data = TRUE)},
  {
    future::plan(sequential)
    set.seed(608, "L'Ecuyer")
    gcipdrtest::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, 
                                        SI_k = 50000, method = 3, checkdata = TRUE,
                                        tabulate.similar.data = TRUE)},
  times = 5
)

print(res)

boxplot(res, names = c("gcipdr", "gcipdrtest"))


### PARALLEL TESTING

# Test 1: Verify parallel reproducibility (same seed, same result across runs)
library(future)
library(future.apply)

# Set parallel backend
future::plan(multisession, workers = 4)

# Run 1
set.seed(608, "L'Ecuyer")
result_par_run1 <- gcipdrtest::Simulate.data.given.IPD(
  testdat, H = 5, stochastic.integration = TRUE,
  SI_k = 50000, method = 3, checkdata = TRUE, tabulate.similar.data = TRUE
)

# Run 2 (same seed)
set.seed(608, "L'Ecuyer")
result_par_run2 <- gcipdrtest::Simulate.data.given.IPD(
  testdat, H = 5, stochastic.integration = TRUE,
  SI_k = 50000, method = 3, checkdata = TRUE, tabulate.similar.data = TRUE
)

# Check: should be identical
all.equal(result_par_run1$Xspace, result_par_run2$Xspace)
# Target: TRUE

# Test 2: Verify sequential still works (backward compatibility)
# Reset to sequential
plan(sequential)

set.seed(608, "L'Ecuyer")
result_seq <- gcipdrtest::Simulate.data.given.IPD(
  testdat, H = 5, stochastic.integration = TRUE,
  SI_k = 50000, method = 3, checkdata = TRUE, tabulate.similar.data = TRUE
)

# Should match the current (non-parallel) benchmark
# Expected time: ~9.48s (same as before)

# Test 3: Benchmark parallel vs sequential


# Sequential
future::plan(sequential)
set.seed(608, "L'Ecuyer")
seq_bench <- microbenchmark::microbenchmark(
  gcipdrtest::Simulate.data.given.IPD(
    testdat, H = 5, stochastic.integration = TRUE,
    SI_k = 50000, method = 3
  ),
  times = 30L
)

# Parallel (4 workers)
future::plan(multisession, workers = 6)
set.seed(608, "L'Ecuyer")
par_bench <- microbenchmark::microbenchmark(
  gcipdrtest::Simulate.data.given.IPD(
    testdat, H = 5, stochastic.integration = TRUE,
    SI_k = 50000, method = 3
  ),
  times = 30L
)

print(seq_bench)
print(par_bench)

# Calculate speedup
speedup <- median(seq_bench$time) / median(par_bench$time)
cat("Parallel speedup (4 workers):", speedup, "x\n")


res <- microbenchmark::microbenchmark(
  {
    future::plan(sequential)
    set.seed(608, "L'Ecuyer")
    gcipdrtest::Simulate.data.given.IPD(
      testdat, H = 5, stochastic.integration = TRUE,
      SI_k = 50000, method = 3
    )
  },
  {
    future::plan(multisession, workers = 4)
    set.seed(608, "L'Ecuyer")
    gcipdrtest::Simulate.data.given.IPD(
      testdat, H = 5, stochastic.integration = TRUE,
      SI_k = 50000, method = 3
    )  
  },
  times = 30L
)


print(res)

boxplot(res, names = c("sequential", "parallel 4"))


# Reset
plan(sequential)