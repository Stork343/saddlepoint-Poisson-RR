##################example1###################
Y.0 = 68
Y.1 = 59
N.0 = 108
N.1 = 124
lam.0 = Y.0/N.0
p.0 = lam.0
lam.1 = Y.1/N.1
p.1 = lam.1
R.hat = lam.0/lam.1
R.hat.var = R.hat^2*( lam.0 / (lam.0 - lam.0^2) + lam.1 / (lam.1 - lam.1^2) ) #计算R估计的方差
R.hat.var = (lam.0/N.0)/lam.1^2 + (lam.0^2 * (lam.1/N.1) ) / lam.1^4
R.hat.var

alpha = 0.05
Z = qnorm(1-0.025,0, 1)
chsiq = qchisq(1 - alpha/2, 1)
cumu3 = function(r){
  n = N.0+N.1
  saddle.func = function(theta){
    lam.0 * exp(theta * N.1) -  lam.1 * r * exp(-theta * N.0 * r)
  }
  theta = uniroot(saddle.func,interval = c(-3, 3), tol=1e-15)$root
  
  phi = exp(N.0 * lam.0 * (exp(theta * N.1) - 1)) *
    exp(N.1 * lam.1 * (exp(-theta * N.0 * r) - 1)) #\varphi_Z(\theta)
  
  mu = N.1 * N.0 * lam.0 * exp(theta * N.1) -
    N.1 * N.0 * lam.1 * r * exp(-theta * N.0 * r) #K_Z'(\theta)
  #phi <- (u / (1 - (1 - u) * exp(theta))) ^ x
  sigma2 = N.1 ^ 2 * N.0  * lam.0 * exp(theta * N.1) +
    N.1 * N.0 ^ 2 * lam.1 * r ^ 2 * exp(-theta * N.0 * r) #K_Z''(\theta)
  
  K3 = N.1 ^ 3 * N.0  * lam.0 * exp(theta * N.1) -
    N.1 * N.0 ^ 3 * lam.1 * r ^ 3 * exp(-theta * N.1 * r) #K_Z^{(3)}(\theta)
  
  K4 = N.1 ^ 4 * N.0  * lam.0 * exp(theta * N.1) +
    N.1 * N.0 ^ 4 * lam.1 * r ^ 4 * exp(-theta * N.0 * r) #K_Z^{(4)}(\theta)
  
  Lam = 1 ^ (1 / 2) * abs(theta) * sigma2 ^ (1 / 2) #\lambda
  gam = exp(-abs(theta)) / (1 - exp(-abs(theta))) #\gamma_{\theta}
  
  B0 = Lam * exp(Lam ^ 2 / 2) * (1 - pnorm(Lam))
  if (is.nan(B0)){
    n = 1
    Lam = n ^ (1 / 2) * abs(theta) * sigma2 ^ (1 / 2) #\lambda
    gam = exp(-abs(theta)) / (1 - exp(-abs(theta))) #\gamma_{\theta}
    
    B0 <- Lam * exp(Lam ^ 2 / 2) * (1 - pnorm(Lam))
  }
  B1 = -Lam * (B0 - (2 * pi) ^ (-1 / 2))
  B2 = Lam ^ 2 * (B0 - (2 * pi) ^ (-1 / 2))
  B3 = -(Lam ^ 3 * B0 - (Lam ^ 3 - Lam) * (2 * pi) ^ (-1 / 2))
  B4 = Lam ^ 4 * B0 - (Lam ^ 4 - Lam ^ 2) * (2 * pi) ^ (-1 / 2)
  B5 = -(Lam ^ 5 * B0 - (Lam ^ 5 - Lam ^ 3 + 3 * Lam) * (2 * pi) ^ (-1 / 2))
  B6 = Lam ^ 6 * B0 - (Lam ^ 6 - Lam ^ 4 + 3 * Lam ^ 2) * (2 * pi) ^ (-1 / 2)
  zeta3 = K3 / (sigma2) ^ (3 / 2)
  zeta4 = K4 / (sigma2) ^ (4 / 2)
  
  a = phi ^ n * exp(theta * r)
  bb = (n * sigma2) ^ (1 / 2) * (1 - exp(-abs(theta)))
  A = a / bb
  I1 = B0
  I2 = 1 / sqrt(n) * (1 / sigma2 ^ (1 / 2) * (1 / abs(theta) - gam) * B1 + zeta3 /
                        6 * sign(theta) * B3)
  p1 = (gam * (1 / 2 - 1 / abs(theta)) + gam ^ 2) * B2 / sigma2
  p2 = 1 / sigma2 ^ (1 / 2) * (1 / abs(theta) - gam) * zeta3 / 6 * sign(theta) * B4
  p3 = zeta4 / 24 * B4 + zeta3 ^ 2 / 72 * B6
  I3 = 1 / n * (p1 + p2 + p3)
  SPA = A * (I1 + I2 + I3)
  return(SPA)
} #SPA main function
f = function(u){
  cumu3(u) - alpha/2
} # find F(R>r)

SPA.L = uniroot(f, interval = c(0.00001, R.hat-0.00001), tol = 1e-15)$root

SPA.U = uniroot(f, interval = c(R.hat+0.00001, 1.5*R.hat),tol = 1e-15)$root

Wald.L = max((R.hat - Z * sqrt(R.hat.var)), 0)
Wald.U = R.hat + Z * sqrt(R.hat.var)


Woolf.L = max(exp(log(R.hat) - Z * (1 / sum(Y.1) + 1 / (N.1 - sum(Y.1)) + 
                                      1 / sum(Y.0) + 1 / (N.0 - sum(Y.0))) ^ 0.5), 0)
Woolf.U = exp(log(R.hat) + Z * (1 / sum(Y.1) + 1 / (N.1 - sum(Y.1)) + 1 /
                                  sum(Y.0) + 1 / (N.0 - sum(Y.0))) ^ 0.5)


LRT.L = max((R.hat - sqrt(chsiq*R.hat.var)), 0)
LRT.U = R.hat + sqrt(chsiq*R.hat.var)  

log.L = max((exp(log(p.0/p.1) - Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + 
                                          ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))), 0)
log.U = min((exp(log(p.0/p.1) + Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + 
                                          ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))), 3*R.hat)

E = lam.1 ^ (2 / 3) - Z ^ 2 * (1 - lam.1) / (9 * N.1 * lam.1 ^ (1 / 3))
D = (lam.1 * lam.0) ^ (1 / 3)
F = lam.0 ^ (2 / 3) - Z ^ 2 * (1 - lam.0) / (9 * N.0 * lam.0 ^ (1 / 3))
F.L = max((((D - sqrt(D ^ 2 - E * F)) / E) ^ 3), 0)
F.U = min((((D + sqrt(D ^ 2 - E * F )) / E) ^ 3), 3 * R.hat)


Exact = poisson.test(c(sum(Y.0),sum(Y.1)), T = c(N.0, N.1), r = R.hat,
                     alternative = "two.sided",
                     conf.level = 0.95)
Exact.L = Exact$conf.int[1]
Exact.U = Exact$conf.int[2]


res = data.frame(
  Low = c(SPA.L, Wald.L, Woolf.L, LRT.L, log.L, F.L, Exact.L),
  Up = c(SPA.U, Wald.U, Woolf.U, LRT.U, log.U, F.U, Exact.U),
  len = c((SPA.U - SPA.L), (Wald.U - Wald.L), 
               (Woolf.U - Woolf.L), (LRT.U - LRT.L),
               (log.U - log.L), (F.U - F.L),
               (Exact.U - Exact.L)),
  Rhat = R.hat )
rownames(res) <- c("SPA", "Wald", 'Woolf', 'LRT', 'log', 'Fieller', "Exact")
res
##################example2###################
Y.0 = 9
Y.1 = 20
N.0 = 50
N.1 = 49
lam.0 = Y.0/N.0
p.0 = lam.0
lam.1 = Y.1/N.1
p.1 = lam.1
R.hat = lam.0/lam.1
R.hat.var = R.hat^2*( lam.0 / (lam.0 - lam.0^2) + lam.1 / (lam.1 - lam.1^2) ) #计算R估计的方差
R.hat.var = (lam.0/N.0)/lam.1^2 + (lam.0^2 * (lam.1/N.1) ) / lam.1^4
R.hat.var

SPA.L = uniroot(f, interval = c(0.00001, R.hat-0.00001), tol = 1e-15)$root

SPA.U = uniroot(f, interval = c(R.hat+0.00001, 1.5*R.hat),tol = 1e-15)$root

Wald.L = max((R.hat - Z * sqrt(R.hat.var)), 0)
Wald.U = R.hat + Z * sqrt(R.hat.var)


Woolf.L = max(exp(log(R.hat) - Z * (1 / sum(Y.1) + 1 / (N.1 - sum(Y.1)) + 
                                      1 / sum(Y.0) + 1 / (N.0 - sum(Y.0))) ^ 0.5), 0)
Woolf.U = exp(log(R.hat) + Z * (1 / sum(Y.1) + 1 / (N.1 - sum(Y.1)) + 1 /
                                  sum(Y.0) + 1 / (N.0 - sum(Y.0))) ^ 0.5)


LRT.L = max((R.hat - sqrt(chsiq*R.hat.var)), 0)
LRT.U = R.hat + sqrt(chsiq*R.hat.var)  

log.L = max((exp(log(p.0/p.1) - Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + 
                                          ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))), 0)
log.U = min((exp(log(p.0/p.1) + Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + 
                                          ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))), 3*R.hat)

E = lam.1 ^ (2 / 3) - Z ^ 2 * (1 - lam.1) / (9 * N.1 * lam.1 ^ (1 / 3))
D = (lam.1 * lam.0) ^ (1 / 3)
F = lam.0 ^ (2 / 3) - Z ^ 2 * (1 - lam.0) / (9 * N.0 * lam.0 ^ (1 / 3))
F.L = max((((D - sqrt(D ^ 2 - E * F)) / E) ^ 3), 0)
F.U = min((((D + sqrt(D ^ 2 - E * F )) / E) ^ 3), 3 * R.hat)


Exact = poisson.test(c(sum(Y.0),sum(Y.1)), T = c(N.0, N.1), r = R.hat,
                     alternative = "two.sided",
                     conf.level = 0.95)
Exact.L = Exact$conf.int[1]
Exact.U = Exact$conf.int[2]


res = data.frame(
  Low = c(SPA.L, Wald.L, Woolf.L, LRT.L, log.L, F.L, Exact.L),
  Up = c(SPA.U, Wald.U, Woolf.U, LRT.U, log.U, F.U, Exact.U),
  len = c((SPA.U - SPA.L), (Wald.U - Wald.L), 
          (Woolf.U - Woolf.L), (LRT.U - LRT.L),
          (log.U - log.L), (F.U - F.L),
          (Exact.U - Exact.L)),
  Rhat = R.hat )
rownames(res) <- c("SPA", "Wald", 'Woolf', 'LRT', 'log', 'Fieller', "Exact")
res


##################example3###################
Y.0 = 83
Y.1 = 3
N.0 = 155
N.1 = 17
lam.0 = Y.0/N.0
p.0 = lam.0
lam.1 = Y.1/N.1
p.1 = lam.1
R.hat = lam.0/lam.1
R.hat.var = R.hat^2*( lam.0 / (lam.0 - lam.0^2) + lam.1 / (lam.1 - lam.1^2) ) #计算R估计的方差
R.hat.var = (lam.0/N.0)/lam.1^2 + (lam.0^2 * (lam.1/N.1) ) / lam.1^4
R.hat.var

SPA.L = uniroot(f, interval = c(0.00001, R.hat-0.00001), tol = 1e-15)$root

SPA.U = uniroot(f, interval = c(R.hat+0.00001, 1.5*R.hat),tol = 1e-15)$root

Wald.L = max((R.hat - Z * sqrt(R.hat.var)), 0)
Wald.U = R.hat + Z * sqrt(R.hat.var)


Woolf.L = max(exp(log(R.hat) - Z * (1 / sum(Y.1) + 1 / (N.1 - sum(Y.1)) + 
                                      1 / sum(Y.0) + 1 / (N.0 - sum(Y.0))) ^ 0.5), 0)
Woolf.U = exp(log(R.hat) + Z * (1 / sum(Y.1) + 1 / (N.1 - sum(Y.1)) + 1 /
                                  sum(Y.0) + 1 / (N.0 - sum(Y.0))) ^ 0.5)


LRT.L = max((R.hat - sqrt(chsiq*R.hat.var)), 0)
LRT.U = R.hat + sqrt(chsiq*R.hat.var)  

log.L = max((exp(log(p.0/p.1) - Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + 
                                          ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))), 0)
log.U = min((exp(log(p.0/p.1) + Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + 
                                          ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))), 3*R.hat)

E = lam.1 ^ (2 / 3) - Z ^ 2 * (1 - lam.1) / (9 * N.1 * lam.1 ^ (1 / 3))
D = (lam.1 * lam.0) ^ (1 / 3)
F = lam.0 ^ (2 / 3) - Z ^ 2 * (1 - lam.0) / (9 * N.0 * lam.0 ^ (1 / 3))
F.L = max((((D - sqrt(D ^ 2 - E * F)) / E) ^ 3), 0)
F.U = min((((D + sqrt(D ^ 2 - E * F )) / E) ^ 3), 3 * R.hat)


Exact = poisson.test(c(sum(Y.0),sum(Y.1)), T = c(N.0, N.1), r = R.hat,
                     alternative = "two.sided",
                     conf.level = 0.95)
Exact.L = Exact$conf.int[1]
Exact.U = Exact$conf.int[2]


res = data.frame(
  Low = c(SPA.L, Wald.L, Woolf.L, LRT.L, log.L, F.L, Exact.L),
  Up = c(SPA.U, Wald.U, Woolf.U, LRT.U, log.U, F.U, Exact.U),
  len = c((SPA.U - SPA.L), (Wald.U - Wald.L), 
          (Woolf.U - Woolf.L), (LRT.U - LRT.L),
          (log.U - log.L), (F.U - F.L),
          (Exact.U - Exact.L)),
  Rhat = R.hat )
rownames(res) <- c("SPA", "Wald", 'Woolf', 'LRT', 'log', 'Fieller', "Exact")
res































