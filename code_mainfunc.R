# 常数值
lam.0 <- 1
lam.1 <- 2
N.0 <- 100
N.1 <- 150

# 定义theta函数
theta <- function(r) {
  (log(r) + log(lam.0/lam.1)) / (N.0 + N.0*r + N.1)
}

# 绘制图像
curve(theta, from = 0.001, to = 2, col = "blue", xlab = "r", ylab = "theta", main = "Graph of theta(r)")

library(rootSolve)
library(extraDistr)

#rtpois(n=10, lambda=0.01, a = 0, b = Inf)



N.1 = 288
N.0 = 289
N.1 = 10
N.0 = 10
lambda.0 = 101/N.1
lambda.1 = 23/N.0

N.1 = 52
N.0 = 36
lam.0 = 23/N.1
lam.1 = 10/N.0
lambda.0 = 0.5
lambda.1 = 0.1
lam.0 = lambda.0
lam.1 = lambda.1
Y.0 = rpois(n=10,lambda = 5)
Y.1 = rpois(n=N.1,lambda = lambda.1)
Y.0
Y.1
lam.0 = mean(Y.0)
lam.0
lam.1 = mean(Y.1)
lam.0/lam.1

loglik <- function(lambda, data) {
  sum(dpois(data, lambda, log = TRUE))
}

library(maxLik)
lam.0 = maxLik(loglik, start = 0, data = Y.0)$estimate
lam.1 = maxLik(loglik, start = 1, data = Y.1)$estimate
R = lam.0/lam.1
R
lam.0.var = 1/mean(Y.0/lam.0^2)
lam.0.var
lam.1.var = 1/mean(Y.1/lam.1^2)
lam.1.var
#0, 0.4662037
n = N.0+N.1
n = 1

-lam.0/lam.1*log(lam.1/lam.0)

for (i in seq(0,2, 0.01)){
  r = i
  saddle.func = function(theta){
    r = i
    N.1 * N.0 * lam.1 * exp(theta * N.0) - N.1 * N.0 * lam.1 * exp(-theta * N.0 * r)
    }
  theta = uniroot(saddle.func,interval = c(-2,2),tol=1e-10)$root
  #theta = N.1 * N.0 * lam.1 * exp(theta * N.0) - N.1 * N.0 * lam.1 * exp(-theta * N.0 * r)
  phi <- exp(N.1*lam.1*(exp(theta*N.0)-1))*
         exp(N.1*lam.1*(exp(-theta*N.0*r)-1)) #\varphi_Z(\theta)
  
  mu <- N.1 * N.0 * lam.1 * exp(theta * N.0) - 
        N.1 * N.0 * lam.1 * r * exp(-theta * N.0 * r) #K_Z'(\theta)
  #phi <- (u / (1 - (1 - u) * exp(theta))) ^ x
  sigma2 <- N.1 * N.0 ^ 2 * lam.1 * exp(theta * N.0) - 
            N.1 * N.0 ^ 2 * lam.1 * r ^ 2 * exp(-theta * N.0 * r) #K_Z''(\theta)
  
  K3 <- N.1 * N.0 ^ 3 * lam.1 * exp(theta * N.0) +
        N.1 * N.0 ^ 3 * lam.1 * r ^ 3 * exp(-theta * N.0 * r) #K_Z^{(3)}(\theta)
  
  K4 <- N.1 * N.0 ^ 4 * lam.1 * exp(theta * N.0) +
        N.1 * N.0 ^ 4 * lam.1 * r ^ 4 * exp(-theta * N.0 * r) #K_Z^{(4)}(\theta)
  
  Lam <- n ^ (1 / 2) * abs(theta) * sigma2 ^ (1 / 2) #\lambda
  gam <- exp(-abs(theta)) / (1 - exp(-abs(theta))) #\gamma_{\theta}
  
  B0 <- Lam * exp(Lam ^ 2 / 2) * (1 - pnorm(Lam))
  B1 <- -Lam * (B0 - (2 * pi) ^ (-1 / 2))
  B2 <- Lam ^ 2 * (B0 - (2 * pi) ^ (-1 / 2))
  B3 <- -(Lam ^ 3 * B0 - (Lam ^ 3 - Lam) * (2 * pi) ^ (-1 / 2))
  B4 <- Lam ^ 4 * B0 - (Lam ^ 4 - Lam ^ 2) * (2 * pi) ^ (-1 / 2)
  B5 <- -(Lam ^ 5 * B0 - (Lam ^ 5 - Lam ^ 3 + 3 * Lam) * (2 * pi) ^ (-1 / 2))
  B6 <- Lam ^ 6 * B0 - (Lam ^ 6 - Lam ^ 4 + 3 * Lam ^ 2) * (2 * pi) ^ (-1 / 2)
  
  zeta3 <- K3 / (sigma2) ^ (3 / 2)
  zeta4 <- K4 / (sigma2) ^ (4 / 2)
  
  
  a <- phi ^ n * exp(-n * theta * r)
  bb <- (n * sigma2) ^ (1 / 2) * (1 - exp(-abs(theta)))
  A <- a / bb
  I1 <- B0
  I2 <- 1 / sqrt(n) * (1 / sigma2 ^ (1 / 2) * (1 / abs(theta) - gam) * B1 + zeta3 /
                         6 * sign(theta) * B3)
  p1 <- (gam * (1 / 2 - 1 / abs(theta)) + gam ^ 2) * B2 / sigma2
  p2 <-
    1 / sigma2 ^ (1 / 2) * (1 / abs(theta) - gam) * zeta3 / 6 *
    sign(theta) * B4
  p3 <- zeta4 / 24 * B4 + zeta3 ^ 2 / 72 * B6
  I3 <- 1 / n * (p1 + p2 + p3)
  sdp <- A * (I1 + I2 + I3)
  #return(sdp)
  print(cat("n:",n,"theta:",theta, "sdp:",sdp,"r",r))
}

cumu = function(r){
#  saddle.func = function(theta){
#    r = r
#    N.1 * N.0 * lam.1 * exp(theta * N.0) - N.1 * N.0 * lam.1 * r * exp(-theta * N.0 * r)
    #N.1 * N.0 * lam.1 * exp(theta * N.0) - N.1 * N.0 * theta * r
#  }
#  theta = uniroot.all(saddle.func,interval = c(-1,1),tol=1e-10)
  #theta= -lam.0/lam.1*log(lam.1/lam.0)
#theta = (log(r))/(N.0+N.0*r)
theta = (log(r)+log(lam.0/lam.1))/(N.0+N.0*r+N.1)
phi <- exp(N.1*lam.1*(exp(theta*N.0)-1))*
  exp(N.1*lam.1*(exp(-theta*N.0*r)-1)) #\varphi_Z(\theta)

mu <- N.1 * N.0 * lam.1 * exp(theta * N.0) - 
  N.1 * N.0 * lam.1 * r * exp(-theta * N.0 * r) #K_Z'(\theta)
#phi <- (u / (1 - (1 - u) * exp(theta))) ^ x
sigma2 <- N.1 * N.0 ^ 2 * lam.1 * exp(theta * N.0) + 
          N.1 * N.0 ^ 2 * lam.1 * r ^ 2 * exp(-theta * N.0 * r) #K_Z''(\theta)

K3 <- N.1 * N.0 ^ 3 * lam.1 * exp(theta * N.0) -
      N.1 * N.0 ^ 3 * lam.1 * r ^ 3 * exp(-theta * N.0 * r) #K_Z^{(3)}(\theta)

K4 <- N.1 * N.0 ^ 4 * lam.1 * exp(theta * N.0) +
      N.1 * N.0 ^ 4 * lam.1 * r ^ 4 * exp(-theta * N.0 * r) #K_Z^{(4)}(\theta)

Lam <- n ^ (1 / 2) * abs(theta) * sigma2 ^ (1 / 2) #\lambda
gam <- exp(-abs(theta)) / (1 - exp(-abs(theta))) #\gamma_{\theta}

B0 <- Lam * exp(Lam ^ 2 / 2) * (1 - pnorm(Lam))
B1 <- -Lam * (B0 - (2 * pi) ^ (-1 / 2))
B2 <- Lam ^ 2 * (B0 - (2 * pi) ^ (-1 / 2))
B3 <- -(Lam ^ 3 * B0 - (Lam ^ 3 - Lam) * (2 * pi) ^ (-1 / 2))
B4 <- Lam ^ 4 * B0 - (Lam ^ 4 - Lam ^ 2) * (2 * pi) ^ (-1 / 2)
B5 <- -(Lam ^ 5 * B0 - (Lam ^ 5 - Lam ^ 3 + 3 * Lam) * (2 * pi) ^ (-1 / 2))
B6 <- Lam ^ 6 * B0 - (Lam ^ 6 - Lam ^ 4 + 3 * Lam ^ 2) * (2 * pi) ^ (-1 / 2)

zeta3 <- K3 / (sigma2) ^ (3 / 2)
zeta4 <- K4 / (sigma2) ^ (4 / 2)


a <- phi^n * exp(theta * r)
bb <- (n * sigma2) ^ (1 / 2) * (1 - exp(-abs(theta)))
A <- a / bb
I1 <- B0
I2 <- 1 / sqrt(n) * (1 / sigma2 ^ (1 / 2) * (1 / abs(theta) - gam) * B1 + zeta3 /
                       6 * sign(theta) * B3)
p1 <- (gam * (1 / 2 - 1 / abs(theta)) + gam ^ 2) * B2 / sigma2
p2 <-
  1 / sigma2 ^ (1 / 2) * (1 / abs(theta) - gam) * zeta3 / 6 *
  sign(theta) * B4
p3 <- zeta4 / 24 * B4 + zeta3 ^ 2 / 72 * B6
I3 <- 1 / n * (p1 + p2 + p3)
sdp <- A * (I1 + I2 + I3)
return(list(sdp=sdp,theta=theta))
}

cumu2 = function(r){
  #  saddle.func = function(theta){
  #    r = r
  #    N.1 * N.0 * lam.1 * exp(theta * N.0) - N.1 * N.0 * lam.1 * r * exp(-theta * N.0 * r)
  #N.1 * N.0 * lam.1 * exp(theta * N.0) - N.1 * N.0 * theta * r
  #  }
  #  theta = uniroot.all(saddle.func,interval = c(-1,1),tol=1e-10)
  #theta= -lam.0/lam.1*log(lam.1/lam.0)
  #theta = (log(r))/(N.0+N.0*r)
  theta = (log(r)+log(lam.0/lam.1))/(N.0+N.0*r+N.1)
  phi <- exp(N.1*lam.1*(exp(theta*N.0)-1))*
    exp(N.1*lam.1*(exp(-theta*N.0*r)-1)) #\varphi_Z(\theta)
  
  mu <- N.1 * N.0 * lam.1 * exp(theta * N.0) - 
    N.1 * N.0 * lam.1 * r * exp(-theta * N.0 * r) #K_Z'(\theta)
  #phi <- (u / (1 - (1 - u) * exp(theta))) ^ x
  sigma2 <- N.1 * N.0 ^ 2 * lam.1 * exp(theta * N.0) + 
    N.1 * N.0 ^ 2 * lam.1 * r ^ 2 * exp(-theta * N.0 * r) #K_Z''(\theta)
  
  K3 <- N.1 * N.0 ^ 3 * lam.1 * exp(theta * N.0) -
    N.1 * N.0 ^ 3 * lam.1 * r ^ 3 * exp(-theta * N.0 * r) #K_Z^{(3)}(\theta)
  
  K4 <- N.1 * N.0 ^ 4 * lam.1 * exp(theta * N.0) +
    N.1 * N.0 ^ 4 * lam.1 * r ^ 4 * exp(-theta * N.0 * r) #K_Z^{(4)}(\theta)
  
  Lam <- n ^ (1 / 2) * abs(theta) * sigma2 ^ (1 / 2) #\lambda
  gam <- exp(-abs(theta)) / (1 - exp(-abs(theta))) #\gamma_{\theta}
  
  B0 <- Lam * exp(Lam ^ 2 / 2) * (1 - pnorm(Lam))
  B1 <- -Lam * (B0 - (2 * pi) ^ (-1 / 2))
  B2 <- Lam ^ 2 * (B0 - (2 * pi) ^ (-1 / 2))
  B3 <- -(Lam ^ 3 * B0 - (Lam ^ 3 - Lam) * (2 * pi) ^ (-1 / 2))
  B4 <- Lam ^ 4 * B0 - (Lam ^ 4 - Lam ^ 2) * (2 * pi) ^ (-1 / 2)
  B5 <- -(Lam ^ 5 * B0 - (Lam ^ 5 - Lam ^ 3 + 3 * Lam) * (2 * pi) ^ (-1 / 2))
  B6 <- Lam ^ 6 * B0 - (Lam ^ 6 - Lam ^ 4 + 3 * Lam ^ 2) * (2 * pi) ^ (-1 / 2)
  
  zeta3 <- K3 / (sigma2) ^ (3 / 2)
  zeta4 <- K4 / (sigma2) ^ (4 / 2)
  
  
  a <- phi^n * exp(theta * r)
  bb <- (n * sigma2) ^ (1 / 2) * (1 - exp(-abs(theta)))
  A <- a / bb
  I1 <- B0
  I2 <- 1 / sqrt(n) * (1 / sigma2 ^ (1 / 2) * (1 / abs(theta) - gam) * B1 + zeta3 /
                         6 * sign(theta) * B3)
  p1 <- (gam * (1 / 2 - 1 / abs(theta)) + gam ^ 2) * B2 / sigma2
  p2 <-
    1 / sigma2 ^ (1 / 2) * (1 / abs(theta) - gam) * zeta3 / 6 *
    sign(theta) * B4
  p3 <- zeta4 / 24 * B4 + zeta3 ^ 2 / 72 * B6
  I3 <- 1 / n * (p1 + p2 + p3)
  sdp <- A * (I1 + I2 + I3)
  return(sdp)
}
#n = N.1
n = 20
for (i in seq(0.001,1,0.001)){
  res = cumu(r=i)$sdp-(alpha/2)
  if (any(is.na(res))) {
    print("存在缺失值")
  } else {
    # 在没有缺失值的情况下进行条件判断
    if (any(abs(res) < 0.001)) {
      #print("存在小于0.001的值")
      #cat("res:",res,"i:",i,"saddlepoint:",cumu(i)$theta,'\n')
      #break
    }else {
    }
  }
  cat("res:",res,"i:",i,"saddlepoint:",cumu(i)$theta,'\n')
}

for (i in seq(0,1,0.001)){
  res = cumu(r=i)$sdp-(alpha/2)
  if (any(is.na(res))) {
    print("存在缺失值")
  } else {
    # 在没有缺失值的情况下进行条件判断
    if (any(abs(res) < 0.001)) {
      print("存在小于0.001的值")
      cat("res:",res,"i:",i,"saddlepoint:",cumu(i)$theta,'\n')
      break
    }else {}
  }
  cat("res:",res,"i:",i,"saddlepoint:",cumu(i)$theta,'\n')
}

f<-function(u)
{cumu2(u)-alpha/2}

g<-function(u)
{cumu2(u)-(1-alpha/2)}

uniroot.all(f, interval = c(0.001,10),tol=1e-10)

uniroot.all(g, interval = c(1,3),tol=1e-3)

R=lam.0/lam.1
R


#wald
#lam.0.var = lam.0/N.0
#lam.1.var = lam.1/N.1
alpha = 0.05
Z = qnorm(1 - alpha/2)
R.var = R^2*( lam.0/(lam.0-lam.0^2) + lam.1/(lam.1-lam.1^2) ) 
R.var = R^2*( lam.0.var/(lam.0.var-lam.0^2) + lam.1.var/(lam.1.var-lam.1^2))

R.var

wald.L = R - Z * sqrt(R.var)
wald.U = R + Z * sqrt(R.var)
wald.L
wald.U
wald.U - wald.L


#log
p.0 = sum(Y.0)/N.0
p.0 = lam.0
p.1 = sum(Y.1)/N.1
p.1 = lam.1
p.0/p.1 #R
log.L = exp(log(p.0/p.1) - Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))
log.U = exp(log(p.0/p.1) + Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))
print(cat("low:",log.L,'up:',log.U ,'length:',(log.U - log.L)))

#待考证
log.lam.0 = log(lam.0)
log.lam.1 = log(lam.1)
log.R = log(lam.0)/log(lam.1)
log.lam.0.var = (1/lam.0)^2*lam.0.var
log.lam.1.var = (1/lam.1)^2*lam.1.var

log.R.var = log.R^2*( log.lam.0.var/(log.lam.0.var-log.lam.0^2) + log.lam.1.var/(log.lam.1.var-log.lam.1^2))
log.R.var
exp(log.R - Z*sqrt(log.R.var))
exp(log.R + Z*sqrt(log.R.var/(N.1+N.0)))



# Woolf
Woolf.L = exp(log(R) - Z*(1/sum(Y.1)+1/(N.1-sum(Y.1))+ 1/sum(Y.0)+1/(N.0-sum(Y.0)) )^0.5)
Woolf.U = exp(log(R) + Z*(1/sum(Y.1)+1/(N.1-sum(Y.1))+ 1/sum(Y.0)+1/(N.0-sum(Y.0)) )^0.5)
print(cat("low:",Woolf.L,'up:',Woolf.U,'length:',(Woolf.U-Woolf.L)))


#LRT
chsiq =  qchisq(1-alpha/2, 1)
LRT.L = R - sqrt(chsiq*R.var)
LRT.U = R + sqrt(chsiq*R.var)
print(cat("low:", LRT.L,'up:', LRT.U, 'length:', (LRT.U-LRT.L)))


#Fieller定理
E = lam.1^(2/3) - Z^2*(1-lam.1)/(9*N.1*lam.1^(1/3))
D = (lam.1*lam.0)^(1/3)
F = lam.0^(2/3) - Z^2*(1-lam.0)/(9*N.0*lam.0^(1/3))
F.L = ((D-sqrt(D^2-E*F))/E)^3
F.U = ((D+sqrt(D^2-E*F))/E)^3
print(cat("low:",F.L,'up:',F.U,'length:',(F.U-F.L)))

