library(sna)
library(BSPADATA)
library(dplyr)
library(ggmcmc)


set.seed(29486173)



binmat.f <- function(n) {
  
  B <- matrix(0, n, n)
  
  B[1, 2] <- 1
  B[n, n-1] <- 1
  
  for (i in 2:(n-1)) {
    
    B[i, i+1] <- 1
    B[i, i-1] <- 1
    
  }
  
  return(B)
  
}

rmvnorm.f <- function(row, mu, Sigma) {
  
  p <- length(mu)
  Q <- chol(Sigma)
  
  Z <- matrix(rnorm(row * p), row, p)
  X <- Z %*% Q + rep(1, row) %*% t(mu)
  X <- data.frame(X)
  
  return(X)
  
}

dmvnorm.f <- function(x, mu, Sigma, log = TRUE) {
  
  p <- length(x)
  
  if (log == !TRUE) {
    
    f <- exp(-0.5 * t(x - mu) %*% solve(Sigma) %*% (x - mu))
    f <- f / ((2 * pi)^(0.5 * p) * det(Sigma)^0.5)
    
  } else {
    
    f <- -0.5 * (log(det(Sigma)) + p * log(2 * pi))
    f <- f - 0.5 * (t(x - mu) %*% solve(Sigma) %*% (x - mu))
    
  }
  
  return(f)
  
}

# invIrW.f <- function() {
#   
#   Function for computing simultaneous autoregressive generating operators (A & A_k)
# 
# }

mu.f <- function(beta, sigma2, rho, omit = TRUE) {
  
  if (omit == !TRUE) {
    
    A <- diag(n) - rho * W
    A <- solve(A)
    
    mu <- A %*% X %*% beta
    
  } else {
    
    A_k <- diag(n-1)
    A_k <- A_k - rho^2 * D_k %*% W %*% s_k %*% t(s_k) %*% W %*% t(D_k) - rho * D_k %*% W %*% t(D_k)
    A_k <- solve(A_k)
    
    X_k <- rho * D_k %*% W %*% s_k %*% t(s_k) %*% X
    X_k <- X_k + D_k %*% X
    
    mu <- A_k %*% X_k %*% beta
    
  }
  
  return(mu)
  
}

Sigma.f <- function(beta, sigma2, rho, omit = TRUE) {
  
  if (omit == !TRUE) {
    
    A <- diag(n) - rho * W
    A <- solve(A)
    
    V <- sigma2 * diag(n)
    
    Sigma <- A %*% V %*% t(A)
    
  } else {
    
    A_k <- diag(n-1)
    A_k <- A_k - rho^2 * D_k %*% W %*% s_k %*% t(s_k) %*% W %*% t(D_k) - rho * D_k %*% W %*% t(D_k)
    A_k <- solve(A_k)
    
    V_k <- rho^2 * D_k %*% W %*% s_k %*% t(s_k) %*% s_k %*% t(s_k) %*% t(W) %*% t(D_k)
    V_k <- V_k + D_k %*% t(D_k)
    V_k <- sigma2 * V_k
    
    Sigma <- A_k %*% V_k %*% t(A_k)
    
  }
  
  return(Sigma)
  
}

delmat.f <- function(k) {
  
  D <- diag(n-1)
  
  D <- as.data.frame(D)
  D <- mutate(D, 0, .before = k)
  D <- as.matrix(D)
  
  return(D)
  
}

selvec.f <- function(k) {
  
  s <- numeric(n)
  
  s[k] <- 1
  
  return(s)
  
}







n <- 50
B <- binmat.f(n = n)

otlr <- 23
# B[(otlr-3):(otlr-2), otlr] <- 1
# B[(otlr+2):(otlr+3), otlr] <- 1
# B[otlr, (otlr+2):(otlr+3)] <- 1
# B[otlr, (otlr-3):(otlr-2)] <- 1
W <- make.stochastic(dat = B, mode = "row")

x0 <- rep(1, n)
x1 <- runif(n, 0, 400)
x2 <- runif(n, 10, 23)
X <- cbind(x0, x1, x2)

rho <- 0.8
beta <- c(18, 0.478, -1.3)
sigma2 <- 45

# A <- invIrW.f()
mu <- mu.f(beta = beta, sigma2 = sigma2, rho = rho, omit = FALSE)
Sigma <- Sigma.f(beta = beta, sigma2 = sigma2, rho = rho, omit = FALSE)
y <- t(rmvnorm.f(row = 1, mu = mu, Sigma = Sigma))

# Type 1 contamination
ctmn <- 5 * sqrt(Sigma[otlr, otlr]) * selvec.f(otlr)
y <- y + ctmn





formula <- y ~ x0 + x1 + x2
data <- data.frame(y = y, x0 = x0, x1 = x1, x2 = x2)

nsim <- 10000
burn <- 2000
step <- 5

prior <- list(b_pri = rep(0, 3), B_pri = diag(1000, 3), r_pri = 0.001, lambda_pri = 0.001)
initial <- list(beta_0 = rep(0, 3), sigma2_0 = 90, rho_0 = 0.5)

MCMC <- hom_sar(formula = formula, data = data, W = W,
                 nsim = nsim, burn = burn, step = step,
                 prior = prior, initial = initial, kernel = "normal")

mcmc <- MCMC$chains
mcmc.df <- ggs(mcmc)
ggmcmc(mcmc.df)





W_IS <- matrix(0, nsim, n)

for (k in 1:n) {
  
  D_k <- delmat.f(k = k)
  s_k <- selvec.f(k = k)
  
  y_k <- D_k %*% y
  
  for (r in 1:nsim) {
    
    beta <- mcmc[r, 1:3]
    sigma2 <- mcmc[r, 4]
    rho <- mcmc[r, 5]
    
    mu <- mu.f(beta = beta, sigma2 = sigma2, rho = rho, omit = FALSE)
    Sigma <- Sigma.f(beta = beta, sigma2 = sigma2, rho = rho, omit = FALSE)
    
    mu_k <- mu.f(beta = beta, sigma2 = sigma2, rho = rho, omit = TRUE)
    Sigma_k <- Sigma.f(beta = beta, sigma2 = sigma2, rho = rho, omit = TRUE)
    
    W_IS[r, k] <- dmvnorm.f(x = y_k, mu = mu_k, Sigma = Sigma_k, log = TRUE)
    W_IS[r, k] <- W_IS[r, k] - dmvnorm.f(x = y, mu = mu, Sigma = Sigma, log = TRUE)
    
  }
  
}



PCA <- prcomp(W_IS, center = TRUE)
summary(PCA)

screeplot(PCA, main = "Scree plot", type = "line")
abline(1, 0, lty = 2)

biplot(PCA, cex=c(0.01, 1), xlab = "PC1 (.9363)", ylab = "PC2 (.03454)")
