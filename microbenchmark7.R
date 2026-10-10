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
library(future)
library(future.apply)
pak::pak("bonorico/gcipdr")

library(gcipdr)


#### START TESTS

testdat <- mtcars[, 1:2]
apply(testdat, 2, mean)
apply(testdat, 2, sd)
SuppDists::moments(testdat$mpg)


test <- gcipdrtest::Simulate.data.given.IPD(
  testdat, H=5, stochastic.integration = TRUE, 
  SI_k = 50000, method = 13, # method 13 user-defined marginals 
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
  SI_k = 50000, method = 13, # method 13 user-defined marginals 
  checkdata = TRUE,
  tabulate.similar.data = TRUE,
  user_defined_marginals = list(
    function(x) qJohnson(x, parms = Jparms1),
    function(x) qnorm(x, mean = 20, sd = 6)
    )
  )




# check kruskal check equal: PASS
res <- microbenchmark::microbenchmark(
  {
    set.seed(608, "L'Ecuyer")
    gcipdr::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, 
                                    SI_k = 50000, method = 7, 
                                    checkdata = TRUE, tabulate.similar.data = TRUE)},
  {
    set.seed(608, "L'Ecuyer")
    gcipdrtest::Simulate.data.given.IPD(testdat, H=5, stochastic.integration = TRUE, 
                                        SI_k = 50000, method = 13, kruskal_use = TRUE,
                                        checkdata = TRUE, tabulate.similar.data = TRUE,
                                        user_defined_marginals = list(
                                          # NOT PASS IF MOMENTS NOT EQUAL TO DECIMAL PLACES
                                          # function(x) gcipdrtest::qgmom(x, mean = 20, sd = 6),
                                          # function(x) gcipdrtest::qgmom(x, mean = 6, sd = 1.8)
                                          function(x) gcipdrtest::qgmom(
                                            x, mean = mean(testdat$mpg), sd = sd(testdat$mpg)),
                                          function(x) gcipdrtest::qgmom(
                                            x, mean = mean(testdat$cyl), sd = sd(testdat$cyl))
                                        )
                                        )
  },
  times = 5,
  check = \(x) {
    
    A <- x[[1]]$similar.data |> 
      map(\(df) df |> 
            as.data.frame()) |> 
      bind_rows()
    B <- x[[2]]$similar.data |> 
      map(\(df) df |> 
            as.data.frame()) |> 
      bind_rows()
    
    all(
      as.vector(as.matrix(A) - as.matrix(B)) |> 
        round(10) == 0
    )
  }
)


## check normal correlation

res <- microbenchmark::microbenchmark(
  {
    set.seed(608, "L'Ecuyer")
    gcipdr::Simulate.data.given.IPD(
      testdat |> 
        mutate(
          bin1 = rbinom(dim(testdat)[1], 1, 0.5)
        ), 
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 3, 
      checkdata = TRUE, tabulate.similar.data = TRUE)},
  {
    set.seed(608, "L'Ecuyer")
    gcipdrtest::Simulate.data.given.IPD(
      testdat |> 
        mutate(
          bin1 = rbinom(dim(testdat)[1], 1, 0.5)
        ), 
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 11, 
      checkdata = TRUE, tabulate.similar.data = TRUE,
      user_defined_marginals = list(
        function(x) gcipdrtest::qgmom(x, mean = 20, sd = 6),
        function(x) gcipdrtest::qgmom(x, mean = 6, sd = 1.8),
        function(x) qbinom(x, 1, 0.5)
      ))
  },
  times = 1
)

#### CHECK EQUAL WITH BINARY MIX

res <- microbenchmark::microbenchmark(
  {
    set.seed(608, "L'Ecuyer")
    gcipdr::Simulate.data.given.IPD(
      testdat |> 
        mutate(
          bin1 = rbinom(dim(testdat)[1], 1, 0.5)
        ), 
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 3, 
      checkdata = TRUE, tabulate.similar.data = TRUE)},
  {
    set.seed(608, "L'Ecuyer")
    dat <- testdat |> 
      mutate(
        bin1 = rbinom(dim(testdat)[1], 1, 0.5)
      )
    gcipdrtest::Simulate.data.given.IPD(
      dat,
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 13, 
      checkdata = TRUE, tabulate.similar.data = TRUE,
      #kruskal_use = TRUE,
      user_defined_marginals = list(
        function(x) gcipdrtest::qgmom(x, mean = mean(dat$mpg), sd = sd(dat$mpg)),
        function(x) gcipdrtest::qgmom(x, mean = mean(dat$cyl), sd = sd(dat$cyl)),
        function(x) qbinom(x, 1, mean(dat$bin1))
      ))
  },
  times = 2,
  check = \(x) {
    
    A <- x[[1]]$similar.data |> 
      map(\(df) df |> 
            as.data.frame()) |> 
      bind_rows()
    B <- x[[2]]$similar.data |> 
      map(\(df) df |> 
            as.data.frame()) |> 
      bind_rows()
    
    all(
      as.vector(as.matrix(A) - as.matrix(B)) |> 
        round(10) == 0
    )
  }
)

## Kruskal USE with binary mix

res <- microbenchmark::microbenchmark(
  {
    set.seed(608, "L'Ecuyer")
    gcipdr::Simulate.data.given.IPD(
      testdat |> 
        mutate(
          bin1 = rbinom(dim(testdat)[1], 1, 0.5)
        ), 
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 3, 
      checkdata = TRUE, tabulate.similar.data = TRUE)},
  {
    set.seed(608, "L'Ecuyer")
    dat <- testdat |> 
      mutate(
        bin1 = rbinom(dim(testdat)[1], 1, 0.5)
      )
    gcipdrtest::Simulate.data.given.IPD(
      dat,
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 13, 
      checkdata = TRUE, tabulate.similar.data = TRUE,
      kruskal_use = TRUE,
      user_defined_marginals = list(
        function(x) gcipdrtest::qgmom(x, mean = mean(dat$mpg), sd = sd(dat$mpg)),
        function(x) gcipdrtest::qgmom(x, mean = mean(dat$cyl), sd = sd(dat$cyl)),
        function(x) qbinom(x, 1, mean(dat$bin1))
      ))
  },
  times = 10
)

boxplot(res, names = c("gcipdr", "gcipdrtest"))


# check parallel


res <- microbenchmark::microbenchmark(
  {
    set.seed(608, "L'Ecuyer")
    gcipdr::Simulate.data.given.IPD(
      testdat |> 
        mutate(
          bin1 = rbinom(dim(testdat)[1], 1, 0.5)
        ), 
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 3, 
      checkdata = TRUE, tabulate.similar.data = TRUE)},
  {
    set.seed(608, "L'Ecuyer")
    dat <- testdat |> 
      mutate(
        bin1 = rbinom(dim(testdat)[1], 1, 0.5)
      )
    gcipdrtest::Simulate.data.given.IPD(
      dat,
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 13, 
      checkdata = TRUE, tabulate.similar.data = TRUE,
      newtrap_parallel = TRUE,
      user_defined_marginals = list(
        function(x) gcipdrtest::qgmom(x, mean = mean(dat$mpg), sd = sd(dat$mpg)),
        function(x) gcipdrtest::qgmom(x, mean = mean(dat$cyl), sd = sd(dat$cyl)),
        function(x) qbinom(x, 1, mean(dat$bin1))
      ))
  },
  times = 10
)

boxplot(res, names = c("gcipdr", "gcipdrtest"))

### add binary var

res <- microbenchmark::microbenchmark(
  {
    set.seed(608, "L'Ecuyer")
    gcipdr::Simulate.data.given.IPD(
      testdat |> 
        mutate(
          bin1 = rbinom(dim(testdat)[1], 1, 0.5),
          bin2 = rbinom(dim(testdat)[1], 1, 0.8)
        ), 
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 3, 
      checkdata = TRUE, tabulate.similar.data = TRUE)},
  {
    set.seed(608, "L'Ecuyer")
    dat <- testdat |> 
      mutate(
        bin1 = rbinom(dim(testdat)[1], 1, 0.5),
        bin2 = rbinom(dim(testdat)[1], 1, 0.8)
      )
    gcipdrtest::Simulate.data.given.IPD(
      dat,
      H=5, stochastic.integration = TRUE, 
      SI_k = 50000, method = 13, 
      checkdata = TRUE, tabulate.similar.data = TRUE,
      newtrap_parallel = TRUE,
      user_defined_marginals = list(
        function(x) gcipdrtest::qgmom(x, mean = mean(dat$mpg), sd = sd(dat$mpg)),
        function(x) gcipdrtest::qgmom(x, mean = mean(dat$cyl), sd = sd(dat$cyl)),
        function(x) qbinom(x, 1, mean(dat$bin1)),
        function(x) qbinom(x, 1, mean(dat$bin2))
      ))
  },
  times = 5
)

boxplot(res, names = c("gcipdr", "gcipdrtest"))
