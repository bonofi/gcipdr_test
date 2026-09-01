# proof of concept test


func_old <- function(H, n) 
  lapply(1:H, function(i)
    rmvnorm(n, mean = rep(0, 5))
  )

func_new <- function(H, n)
  rmvnorm(n*H, mean = rep(0, 5)) |>
  as.data.frame() |> 
  split(
    rep(1:H, rep(n, H))
  )
    
func_new2 <- function(H, n)
  cbind( 
    rmvnorm(n*H, mean = rep(0, 5)),
    did = rep(1:H, rep(n, H))
    )
    



H <- 300
n <- 100
res <- microbenchmark(
  func_old(H,n),
  func_new(H,n),
  func_new2(H,n),
  times = 100L
)


print(res)
boxplot(res, names = c("lapply", "vec+split", "vector."))