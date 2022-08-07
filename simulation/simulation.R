#                                             CL
#                                             29/05/22


#                                             Packages



# install.packages("pracma")
# install.packages("mvtnorm")
# install.packages("ggplot2")
# install.packages("MCMCpack")
# install.packages("ggmcmc")
# install.packages("dplyr")
# install.packages("scatterplot3d")


library(pracma)
library(mvtnorm)
library(ggplot2)
library(MCMCpack)
library(ggmcmc)
library(dplyr)
library(scatterplot3d)



#                                             Replicas


# set.seed(743435)
# set.seed(29486173)

set.seed(80357)
# set.seed(80433)


#                                             Palette


hex = c("#000000", "#FFFFFF", "#FF0000", "#00FF00",
        "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF",
        "#C0C0C0", "#808080", "#800000", "#808000",
        "#008000", "#800080", "#008080", "#000080")



hex_50 = c(rep(hex[4], 20), rep(hex[8], 5), rep(hex[4], 25))
id_50 = c(rep("", 19), paste0("", 20:26), rep("", 24))




#                                             Datum



rho_0 = 0.9

m = 5
zeros_m = numeric(m)      # 0
eye_m = diag(m)           # 1
beta_0 = rep(0.5, m)

sigma2_0 = 0.5




n = 50
zeros_n = numeric(n)      # 0
eye_n = diag(n)           # 1

C = Diag(rep(1, n - 1), 1) + Diag(rep(1, n - 1), -1)
W = C / rowSums(C)      #      right stochastic spatial weight matrix

X = rmvnorm(n, zeros_m, eye_m)    #  matrix normal design matrix

e_0 = rmvnorm(n, numeric(1), sigma2_0 * diag(1))    #  multivariate normal vector of disturbances
A_0 = eye_n - rho_0 * W
y_0 = solve(A_0) %*% X %*% beta_0 + solve(A_0) %*% e_0



#                                             Outliers


function_0 = function(first, last, deviate) {


  # additive outliers

  y_0[seq(first, last)] = y_0[seq(first, last)] + deviate * sqrt(sigma2_0)

  output = y_0

  return(output)


}

first = 21
last = 25

deviate = 5
# deviate = 4.2
# deviate = 0


y = function_0(first, last, deviate)    #   influential observations will distort parameter estimates
plot(y, ylim = c(-5, 5))

y_df = as.data.frame(y)
ggplot(y_df, aes(1:50, V1)) +           #   visual 01
  xlab("Index") +
  ylab("Value") +
  geom_point(color = hex_50, size = 5) +
  theme_bw()


#                                             Estimation




function_03 = function(rho, beta, sigma2) {
  
  # unperturbed loglikelihood
  
  A = eye_n - rho * W
  V = sigma2 * eye_n

  e = A %*% y - X %*% beta

  output = dmvnorm(t(e), zeros_n, V, log = TRUE)  # mvtnorm::dmvnorm vectors in rows
  output = output + log(det(A))
  
  
  
  return(output)
  
}




function_04 = function(rho, beta, sigma2) {
  
  # rho conditional
  

  
  A = eye_n - rho * W
  
  e = A %*% y - X %*% beta
  
  output = exp(-t(e) %*% e / (2 * sigma2))
  output = output * det(A)
  
  return(output)
  
}



function_05 = function(tau, beta, sigma2) {


  # tau conditional

  rho = (exp(tau) * u_0 + l_0) / (1 + exp(tau))
  J = (exp(tau) * (u_0 - l_0)) / (1 + exp(tau))^2
  
  
  output = function_04(rho, beta, sigma2)
  output = output * abs(J)

  return(output)

}



function_06 = function(ratio, candidate, last) {
  
  # accept/ reject
  
  if (runif(1, 0, 1) <= ratio) {
    
    output = candidate
    
  } else {
    
    output = last
    
  }
  
  return(output)
  
}


function_07 = function(rho, beta, sigma2) {


  # rho metropolis

  tau_last = log((rho - l_0) / (u_0 - rho))      #    transform

  tau_candidate = rnorm(1, tau_last, rho_tune)       #     random walk proposal
  logratio = log(function_05(tau_candidate, beta, sigma2))
  logratio = logratio - log(function_05(tau_last, beta, sigma2))
  ratio = min(exp(logratio), 1)
  output = function_06(ratio, tau_candidate, tau_last)

  output = (exp(output) * u_0 + l_0) / (1 + exp(output))   #   back transform

  return(output)

}


# function_07 = function(rho, beta, sigma2) {
# 
#   # rho metropolis (alt.)
# 
#   rho_last = rho
# 
#   rho_candidate = runif(1, l_0, u_0)      #       uniform proposal
# 
#   logratio = log(function_04(rho_candidate, beta, sigma2))
#   logratio = logratio - log(function_04(rho_last, beta, sigma2))
# 
#   ratio = min(exp(logratio), 1)
# 
#   output = function_06(ratio, rho_candidate, rho_last)
# 
#   return(output)
# 
# 
# 
# 
# }


function_08 = function(rho, beta, sigma2) {
  
  # beta gibbs
  
  D_1 = solve(t(X) %*% X + solve(D_0))
  
  A = eye_n - rho * W
  
  
  c_1 = D_1 %*% (t(X) %*% A %*% y + solve(D_0) %*% c_0)
  
  
  output = rmvnorm(1, c_1, sigma2 * D_1)
  
  return(output)
  
}



function_09 = function(rho, beta, sigma2) {
  
  
  # sigma2 gibbs
  
  A = eye_n - rho * W
  e = A %*% y - X %*% beta
  
  b_1 = (t(e) %*% e + t(beta - c_0) %*% solve(D_0) %*% (beta - c_0)) / 2 + b_0
  a_1 = (n + m) / 2 + a_0
  output = rinvgamma(1, a_1, b_1)
  
  return(output)
  
}





l_0 = -1
u_0 = 1
a_0 = 0
b_0 = 0
c_0 = numeric(m)
D_0 = diag(m)



I = 5000  #             iterations
B = 500   #             burn in
R = I - B #             retained draws
t = 1     #             thinning
MCMC = matrix(0, I, m + 2)
MCMC[1, ] = 0.5
rho_tune = 1

for (i in 2:I) {
  
  

  
  
  rho_update = function_07(MCMC[i - 1, 1], MCMC[i - 1, 2:(m + 1)], MCMC[i - 1, m + 2])
  MCMC[i, 1] = rho_update
  

  
  
  beta_update = function_08(MCMC[i, 1], MCMC[i - 1, 2:(m + 1)], MCMC[i - 1, m + 2])
  MCMC[i, 2:(m + 1)] = beta_update
  

  
  sigma2_update = function_09(MCMC[i, 1], MCMC[i, 2:(m + 1)], MCMC[i - 1, m + 2])
  MCMC[i, m + 2] = sigma2_update
  
  
  
  
}





MCMC = MCMC[seq(B + 1, I, t), ]         # discard burn in period & thin markov chain

MCMC = mcmc(MCMC)
colnames(MCMC) = c("rho", paste0("beta", 1:m), "sigma2")
summary(MCMC)
rejectionRate(MCMC)



MCMC_GG = ggs(MCMC)
ggmcmc(MCMC_GG)          #     visual 02





#                                             Case Influence Analysis


D = diag(n - 1)

s = numeric(n)




function_10 = function(k) {

  # deletion matrix


  D = as.data.frame(D)

  output = mutate(D, 0, .before = k)

  output = as.matrix(output)

  return(output)


}



function_11 = function(k) {

  # selection vector

  s[k] = 1
  output = s

  return(output)

}




function_14 = function(rho, beta, sigma2, k) {

  # perturbed loglikelihood


  D_k = function_10(k)
  s_k = function_11(k)

  A_k = diag(n - 1) - rho^2 * D_k %*% W %*% s_k %*% t(s_k) %*% W %*% t(D_k)
  A_k = A_k - rho * D_k %*% W %*% t(D_k)
  V_k = sigma2 * rho^2 * D_k %*% W %*% s_k %*% t(s_k) %*% s_k %*% t(s_k) %*% t(W) %*% t(D_k)
  V_k = V_k + sigma2 * D_k %*% t(D_k)

  y_k = D_k %*% y
  X_k = rho * D_k %*% W %*% s_k %*% t(s_k) %*% X + D_k %*% X
  e_k = A_k %*% y_k - X_k %*% beta


  output = dmvnorm(t(e_k), numeric(n - 1), V_k, log = TRUE)  # mvtnorm::dmvnorm vectors in rows
  output = output + log(det(A_k))



  return(output)

}


function_15 = function(rho, beta, sigma2, k) {

  # unnormalised weight

  output = function_14(rho, beta, sigma2, k)
  output = output - function_03(rho, beta, sigma2)


  return(output)


}


W_IS = matrix(0, R, n)


for (k in 1:n) {



  for (r in 1:R) {



    weight = function_15(MCMC[r, 1], MCMC[r, 2:(m + 1)], MCMC[r, m + 2], k)
    W_IS[r, k] = weight


  }


}

C_IS = cov(W_IS)


#                                             Static Eigen Analysis



PCA = prcomp(C_IS)
PCA_R = PCA$rotation         # rotation matrix


PCA_3D = scatterplot3d(x = PCA_R[, 1], y = PCA_R[, 2], z = PCA_R[, 3],        #  visual 03
                       xlab = "PC1 (.76)", ylab = "PC2 (.13)", zlab = "PC3 (.06)",
                       color = hex_50, pch = 19, type = "h")
PCA_2D = PCA_3D$xyz.convert(PCA_R[,1], PCA_R[,2], PCA_R[,3])
text(PCA_2D$x, PCA_2D$y,
     labels = id_50, pos = 3)



#                                             Interactive Eigen Analysis




