# MICROBENCHMARK TESTS between bonorico/gcipdr and bonofi/gcipdr_test@optimization_Rccp3

custom_check <- function(res, extra_check = FALSE){
  
  tol <- 0.01 # tolerance parameter for correlations
  x <- res[[1]]
  y <- res[[2]]
  
  gc_param_x <- x$copula.parameters
  gc_param_y <- y$copula.parameters
  
  checks <- NA
  if (extra_check)
  {
    # extract differences for each moment 
    checks <- lapply(
      x$is.data.similar |> 
        names() |>
        # except bool
        head(-1),
      function(name)
      {
        # check correlations
        corr_x <- x$is.data.similar[[name]]
        corr_y <- y$is.data.similar[[name]]
        
        corr_diff <- corr_x |> 
          as.data.frame() |> 
          select(diff) |> 
          rename(diff1 = diff) |>
          rownames_to_column() |> 
          inner_join(
            corr_y |> 
              as.data.frame() |> 
              select(diff) |> 
              rename(diff2 = diff) |> 
              rownames_to_column(),
            by = "rowname"
          ) |> rowwise() |> 
          mutate(diff = diff1 - diff2) |> 
          pull(diff)
        
        
        corrnames_diff <- setdiff(
          rownames(corr_x),
          rownames(corr_y)
        )
        if (length(corrnames_diff) > 0)
          warning(
            paste0(
              name, " between ref and new differ by ",
              paste(corrnames_diff, collapse = ", ")
            )
          )
        
        
        message(
          paste0(
            "Mean difference", name, ": ", 
            mean(corr_diff, na.rm = TRUE)
          )
        )  
        
        
        all(abs(corr_diff) <= tol)
        
      }
    )
    
  }
  
  # check copula parameters difference with tolerance value
  diffparam <- (gc_param_x[lower.tri(gc_param_x)] - gc_param_y[lower.tri(gc_param_y)])
  check1 <- all(abs(diffparam) <= tol)
  
  message(
    paste0(
      "Mean difference copula params: ", mean(diffparam, na.rm = TRUE)
    )
  )
  
  all(
    na.omit(
      c(
        check1, 
        unlist(checks),
        x$is.data.similar$bool == y$is.data.similar$bool 
      )
    )
  )
}

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

install_github("bonofi/gcipdr_test", ref="optimization_Rccp3")

library(gcipdr)
library(gcipdrtest)

################### TEST EQUIVALNCE OF C++ function with R for rmvnorm

K   <- 100      # small K for easy visual inspection
rz  <- 0.2
Sigmaz <- matrix(c(1, rz, rz, 1), 2, 2)

# USING CHOLESKY DECOMPOSITION
# Test 1: rmvnorm (original R)
set.seed(42)
z_r <- mvtnorm::rmvnorm(K, sigma = Sigmaz, method = "chol")
print(round(z_r, 6) |> head())

# Test 2: our C++ version
set.seed(42)
z_cpp <- gcipdrtest::test_rmvnorm_cpp_cholesky(K, rz)
print(round(z_cpp, 6) |> head())

# Test 3: are they identical?
cat("Max absolute difference:", max(abs(z_r - z_cpp)), "\n")
# Target: 0 (or machine epsilon ~1e-16)

# USING EIGEN DECOMPOSITION
# Test 1: rmvnorm (original R)
set.seed(42)
z_r_eigen <- mvtnorm::rmvnorm(K, sigma = Sigmaz)
print(round(z_r_eigen, 6) |> head())

# Test 2: our C++ version
set.seed(42)
z_cpp_eigen <- gcipdrtest::test_rmvnorm_cpp_eigen(K, rz)
print(round(z_cpp_eigen, 6) |> head())

# Test 3: are they identical?
cat("Max absolute difference:", max(abs(z_r_eigen - z_cpp_eigen)), "\n")
# Target: 0 (or machine epsilon ~1e-16)



#### START TESTS

testdat <- mtcars[, 1:5]

# Compare if Rccp conversion 
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
    set.seed(608, "L'Ecuyer")
    gcipdrtest::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, 
                                        SI_k = 50000, method = 3, checkdata = TRUE,
                                        tabulate.similar.data = TRUE)},
  times = 30L,
  check = "equal"
)

print(res)

boxplot(res, names = c("gcipdr", "gcipdrtest"))


### ASEESS DOMINANCE OG GX1 GX2

# Test: how much time is spent in Gx1/Gx2 callbacks vs rest?
K <- 50000
rz <- 0.3

# Time the marginal inverse call alone (assuming gamma marginal)
p_test <- runif(K)
system.time({
  for (i in 1:30){
    x1 <- gcipdrtest::qgmom(p_test, 3, 1)
    x2 <- gcipdrtest::qgmom(p_test, 3, 1)
  } 
    
})

# Time the full mccovx1x2_cpp call
system.time({
  for (i in 1:30) mccovx1x2_cpp(
    rz = 0.3, 
    Gx1 = \(x) gcipdrtest::qgmom(x, 3, 1), 
    Gx2 = \(x) gcipdrtest::qgmom(x, 3, 1), 
    rx = 0, 
    meanx = 3, 
    sdx = 1, 
    TRUE, TRUE, K)
})
