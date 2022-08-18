library(sna)
library(RColorBrewer)
library(mvtnorm)
library(spatialreg)




set.seed(743435)




function_01 <- function(n) {
  
  C <- matrix(0, n, n)
  
  C[1, 2] <- 1
  C[n, n-1] <- 1
  
  for (i in 2:(n-1)) {
    
    C[i, i-1] <- 1
    C[i, i+1] <- 1
    
  }
  
  return(C)
  
}

function_02 <- function(data, first, last, deviate) {
  
  y <- data
  
  y[first:last] <- y[first:last] + deviate
  
  return(y)
  
}

# function_mvnd <- function(x, p, mu, Sigma) {
#   
#   f <- exp(-0.5 * t(x - mu) %*% solve(Sigma) %*% (x - mu))
#   f <- f / ((2 * pi)^(p / 2) * det(Sigma)^0.5)
#   
#   return(f)
#   
# }







n <- 50
C <- function_01(n)
W <- sna::make.stochastic(dat = C, mode = "row")

p <- 1
X <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = diag(p))

rho_0 <- 0.5
A_0 <- spatialreg::invIrW(x = W, rho = rho_0)

sigma2_0 <- 0.5
e_0 <- rnorm(n, 0, sqrt(sigma2_0))

beta_0 <- rep(0.5, p)
y_0 <- A_0 %*% X %*% beta_0 + A_0 %*% e_0

y <- function_02(y_0, 25, 25, 3 * sqrt(sigma2_0))





image(1:n, 1:n, t(W[, n:1]),
      breaks = c(0, 0.25, 0.75, 1),
      col = RColorBrewer::brewer.pal(3, "Purples"),
      main = "Spatial weights representation for y",
      axes = FALSE,
      xlab = "", ylab = "")
box()





