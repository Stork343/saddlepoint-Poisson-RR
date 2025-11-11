SPA = function(r) {
  N.0 = N.0
  N.1 = N.1
  lam.0 = lam.0
  lam.1 = lam.1
  r = r
  saddle.func = function(theta) {
    lam.0 * exp(theta * N.1) -  lam.1 * r * exp(-theta * N.0 * r)
  }
  theta = uniroot(saddle.func, interval = c(-5, 5), tol = 1e-15)$root
  mu = N.1 * N.0 * lam.0 * exp(theta * N.1) -
    N.1 * N.0 * lam.1 * r * exp(-theta * N.0 * r) #K_Z'(\theta)
  sigma2 = N.1 ^ 2 * N.0  * lam.0 * exp(theta * N.1) +
    N.1 * N.0 ^ 2 * lam.1 * r ^ 2 * exp(-theta * N.0 * r) #K_Z''(\theta)
  phi = exp(N.0 * lam.0 * (exp(theta * N.1) - 1)) *
    exp(N.1 * lam.1 * (exp(-theta * N.0 * r) - 1))
  w = sign(theta) * sqrt(2 * (theta * r - log(phi)))
  v = theta * sqrt(sigma2)
  #d = (1 / (2 * pi * sigma2)) ^ 0.5 * exp(Kz(N.0, N.1, lam.0, lam.1, r) - theta * r)
  f = pnorm(w) + dnorm(v) * (1 / w - 1 / v)
  return(f)
}
alpha = 0.05
f = function(u) {
  SPA(u) - alpha / 2
}
g = function(u) {
  SPA(u) - (1 - alpha / 2)
}




N.0 = 20
N.1 = 10
lambda.0 = 1
lambda.1 = 3

SPA.L = NULL
SPA.U = NULL

for (i in 1:2000){

Y.0 = rpois(n = N.0, lambda = lambda.0) #sampling
Y.1 = rpois(n = N.1, lambda = lambda.1) #sampling

if (sum(Y.0) >= 1 & sum(Y.1) >= 1){
  lam.0 = mean(Y.0) #MLE
  lam.1 = mean(Y.1) #MLE
}
R.hat = lam.0 / lam.1 #计算R估计
R.hat
cat(R.hat, lam.0, lam.1,'\n')

root1 = uniroot(f, interval = c(0.01, 10), tol = 1e-15)$root
root2 = uniroot(g, interval = c(0.01, 10), tol = 1e-15)$root

root.L = c(root.L, root1)
root.U = c(root.U, root2)
}
mean(root.L)
mean(root.U)


library(MASS)

N.0 = 20
N.1 = 10
lambda.0 = 1
lambda.1 = 3
root.L = NULL 
root.U = NULL
CI = NULL
count = 0
##ratio with CI between two sample means
for (i in 1:2000){
  out <- mvrnorm(10, mu = c(lambda.0, lambda.1), Sigma = matrix(c( lambda.0, 0.2, 0.2, lambda.1),
                                                   ncol = 2))
  res = FiellerRatio(mean(out[, 1]), mean(out[, 2]), V = cov(out)/10)
  root.L = c(root.L, res[2])
  root.U = c(root.U, res[3])
  CI = c(CI, res[3] - res[2])
  if (res[2] < res[1] && res[1] < res[3]){
    count = count + 1
  }
}
mean(root.L)
mean(root.U)
mean(CI)
count/2000


p.0 = lam.0
p.1 = lam.1
log.L = max((exp(log(p.0/p.1) - Z*sqrt( (( abs(N.1 - sum(Y.1)) )/sum(Y.1))/N.1 + 
                                          (( abs(N.0 - sum(Y.0)) )/sum(Y.0))/N.0 ))), 0)
log.U = min((exp(log(p.0/p.1) + Z*sqrt( (( abs(N.1 - sum(Y.1)) )/sum(Y.1))/N.1 + 
                                          (( abs(N.0 - sum(Y.0)) )/sum(Y.0))/N.0 ))), 3*R.hat)
#print(cat("low:",log.L,'up:',log.U ,'length:',(log.U - log.L)))
if (!is.nan(log.L) & !is.nan(log.U) & !is.infinite(log.L) & !is.infinite(log.U)){
  if (log.L <= R & R <= log.U){
    log.K = log.K + 1
    log.Low = c(log.Low, log.L)
    log.Up = c(log.Up, log.U)
  }
}    

log.L
log.U
