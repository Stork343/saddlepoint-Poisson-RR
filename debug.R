alpha=0.05
r = 
cumu2 = function(r) {
  r = r
  saddle.func = function(theta){
    lam.0 * exp(theta * N.1) -  lam.1 * r * exp(-theta * N.0 * r)
  }
  theta = uniroot(saddle.func,interval = c(-10,10), tol=1e-10)$root
  if (is.nan(theta)){
    theta=0.1
  }
  phi <- exp(N.0 * lam.0 * (exp(theta * N.1) - 1)) *
         exp(N.1 * lam.1 * (exp(-theta * N.0 * r) - 1)) #\varphi_Z(\theta)
  
  mu <- N.1 * N.0 * lam.0 * exp(theta * N.1) -
        N.1 * N.0 * lam.1 * r * exp(-theta * N.0 * r) #K_Z'(\theta)
  #phi <- (u / (1 - (1 - u) * exp(theta))) ^ x
  kapp2 <- N.1 ^ 2 * N.0  * lam.0 * exp(theta * N.1) +
            N.1 * N.0 ^ 2 * lam.1 * r ^ 2 * exp(-theta * N.0 * r) #K_Z''(\theta)
  
  rho3 <- N.1 ^ 3 * N.0  * lam.0 * exp(theta * N.1) -
        N.1 * N.0 ^ 3 * lam.1 * r ^ 3 * exp(-theta * N.1 * r) #K_Z^{(3)}(\theta)
  
  K4 <- N.1 ^ 4 * N.0  * lam.0 * exp(theta * N.1) +
        N.1 * N.0 ^ 4 * lam.1 * r ^ 4 * exp(-theta * N.0 * r) #K_Z^{(4)}(\theta)
  
  Lam <- n ^ (1 / 2) * abs(theta) * sigma2 ^ (1 / 2) #\lambda
  gam <- exp(-abs(theta)) / (1 - exp(-abs(theta))) #\gamma_{\theta}
  
  B0 <- Lam * exp(Lam ^ 2 / 2) * (1 - pnorm(Lam))
  B1 <- -Lam * (B0 - (2 * pi) ^ (-1 / 2))
  B2 <- Lam ^ 2 * (B0 - (2 * pi) ^ (-1 / 2))
  B3 <- -(Lam ^ 3 * B0 - (Lam ^ 3 - Lam) * (2 * pi) ^ (-1 / 2))
  B4 <- Lam ^ 4 * B0 - (Lam ^ 4 - Lam ^ 2) * (2 * pi) ^ (-1 / 2)
  B5 <-
    -(Lam ^ 5 * B0 - (Lam ^ 5 - Lam ^ 3 + 3 * Lam) * (2 * pi) ^ (-1 / 2))
  B6 <-
    Lam ^ 6 * B0 - (Lam ^ 6 - Lam ^ 4 + 3 * Lam ^ 2) * (2 * pi) ^ (-1 / 2)
  
  zeta3 <- K3 / (sigma2) ^ (3 / 2)
  zeta4 <- K4 / (sigma2) ^ (4 / 2)
  
  
  a <- phi ^ n * exp(theta * r)
  bb <- (n * sigma2) ^ (1 / 2) * (1 - exp(-abs(theta)))
  A <- a / bb
  I1 <- B0
  I2 <-
    1 / sqrt(n) * (1 / sigma2 ^ (1 / 2) * (1 / abs(theta) - gam) * B1 + zeta3 /
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
f <- function(u){
  return(cumu2(u) - alpha / 2)
} # find F(R>r)


# 二分法求根函数（返回NULL如果无解）
# 二分法求根函数（返回NULL如果无解或函数值不是数值型）
bisection <- function(f, interval, tol = 1e-5, max_iter = 1000) {
  a <- interval[1]
  b <- interval[2]
  
  # 尝试计算初始端点的函数值
  tryCatch({
    f_a <- f(a)
    f_b <- f(b)
  }, error = function(e) {
    cat("Error in evaluating the function at interval endpoints.\n")
    return(NULL)
  })
  
  # 检查函数值是否为数值型
  if (!is.numeric(f_a) || !is.numeric(f_b)) {
    cat("Function values at interval endpoints are not numeric.\n")
    return(NULL)
  }
  
  if (f_a * f_b > 0) {
    stop("The function values at interval endpoints must have opposite signs.")
  }
  
  iter <- 1
  
  while ((b - a) / 2 > tol && iter <= max_iter) {
    c <- (a + b) / 2
    
    # 尝试计算函数值
    tryCatch({
      f_c <- f(c)
    }, error = function(e) {
      cat("Error in evaluating the function at the midpoint.\n")
      return(NULL)
    })
    
    # 检查函数值是否为数值型
    if (!is.numeric(f_c)) {
      cat("Function value at the midpoint is not numeric.\n")
      return(NULL)
    }
    
    if (f_c == 0) {
      return(c)
    } else if (f_c * f_a < 0) {
      b <- c
    } else {
      a <- c
    }
    iter <- iter + 1
  }
  
  # 判断是否找到根
  if (iter > max_iter) {
    warning("Maximum number of iterations reached. No root found.")
    return(NULL)
  }
  
  return((a + b) / 2)
}

lam=1
x = rpois(20, lam)
x
xSeq <- seq(0, 9, 1)
xSeq
tmp <- dsaddle(y = xSeq, X = x, decay = 0.05, log = F,normalize=T)  # Un-normalized EES
tmp2 <- dsaddle(y = xSeq, X = x, decay = 0.05,             # EES normalized by importance sampling
                normalize = TRUE, control = list("method" = "IS", nNorm = 500), log = TRUE)
plot(xSeq, exp(tmp2$llk), type = 'l', ylab = "Density", xlab = "x")
lines(xSeq, dpois(xSeq, lam), col = 3)
#plot(xSeq, dpois(xSeq, lam), type = 'l', ylab = "Density", xlab = "x")

n=100
lam0=0.5
m = 100
lam1=0.1
x = rpois(n,lam0)
y = rpois(m,lam1)
a = poisson.test(c(sum(x),sum(y)), T = c(n,m), r = 5,
                 alternative = "two.sided",
                 conf.level = 0.95)
a$conf.int
a$estimate
#a$method

n=(sum(x)+sum(y))
p = n/(n+m)

b = binom.test(sum(x), n=(sum(x)+sum(y)), n/(n+m), "two.sided")$conf.int
b

b = m*b/(m*(1-b))
b

n = 50
1/(1+z^2/n) * (p+z^2/(2*n)) + z/(1 + z^2/n)*sqrt(p*(1-p)/n + z^2/(4*n^2))

1/(1+z^2/n) * (p+z^2/(2*n)) - z/(1 + z^2/n)*sqrt(p*(1-p)/n + z^2/(4*n^2))

pbinom(sum(x), sum(x)+sum(y), n/(n+m))

qbeta(0.025, sum(x), sum(y)+ 1)
qbeta(1 - alpha, sum(x) + 1, n - x)