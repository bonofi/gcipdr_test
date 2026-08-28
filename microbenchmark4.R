# MICROBENCHMARK TESTS between bonorico/gcipdr and bonofi/gcipdr_test@optimization_Rccp

# install dependencies

url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
pkgFile <- "JohnsonDistribution_0.24.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)

library(tidyverse)
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
out <- c()
for ( i in 1:10){
  
  set.seed(608+i, "L'Ecuyer")
  ref <- gcipdr::Simulate.data.given.IPD(testdat, H=1, stochastic.integration = TRUE, 
                                         SI_k = 50000, method = 3, checkdata = TRUE,
                                         tabulate.similar.data = TRUE)
  
  set.seed(608+i, "L'Ecuyer")
  new <- gcipdrtest::Simulate.data.given.IPD(testdat, H=1, stochastic.integration = TRUE, 
                                             SI_k = 50000, method = 3, checkdata = TRUE,
                                             tabulate.similar.data = TRUE)
  # check copula parameters
  out <- c(
    out,
    custom_check(list(ref, new))
  ) 
}; all(out)


## check extended

### check equivalence with custom check (increase H) on 1-time benchmark
out <- c()
for ( i in 1:10){
  
  set.seed(608+i, "L'Ecuyer")
  ref <- gcipdr::Simulate.data.given.IPD(testdat, H=5000, stochastic.integration = TRUE, 
                                         SI_k = 50000, method = 3, checkdata = TRUE,
                                         tabulate.similar.data = TRUE)
  
  set.seed(608+i, "L'Ecuyer")
  new <- gcipdrtest::Simulate.data.given.IPD(testdat, H=5000, stochastic.integration = TRUE, 
                                             SI_k = 50000, method = 3, checkdata = TRUE,
                                             tabulate.similar.data = TRUE)
  # check copula parameters
  out <- c(
    out,
    custom_check(list(ref, new), extra_check = TRUE)
  ) 
}; all(out)