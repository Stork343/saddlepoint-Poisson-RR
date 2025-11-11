options(scipen = 999, digits = 6) #打印精度设置
library(rootSolve)
library(extraDistr)
library(latex2exp)
simulation = function(sim.N, N.0, N.1, lambda.0, lambda.1, alpha = 0.05) {
# sim.N 模拟次数 
# N.0 暴露组样本量
# N.1 对照组样本量
# lambda.0 暴露组发病率
# lambda.1 对照组发病率
# alpha 置信水平 

#######初始化区间#######
  SPA.Low = NULL
  SPA.Up = NULL
  Wald.Low = NULL
  Wald.Up = NULL
  Woolf.Low = NULL
  Woolf.Up = NULL
  LRT.Low = NULL
  LRT.Up = NULL
  log.Low = NULL
  log.Up = NULL
  F.Low = NULL
  F.Up = NULL
  Norm.Low = NULL
  Norm.Up = NULL
  Exact.Low = NULL
  Exact.Up = NULL
  Rhat = NULL
  SPA.K = 0
  Wald.K = 0
  Woolf.K = 0
  LRT.K = 0
  log.K = 0
  F.K = 0
  Norm.K = 0
  Exact.K = 0
  fault = 0
  root1 = NULL
  root2 = NULL
  #########分位点初始化#########
  alpha = alpha #名义水平
  Z = qnorm(1 - alpha/2) #正态临界
  chsiq = qchisq(1 - alpha/2, 1) #卡方临界
  
  ##########参数初始化##########
  R = lambda.0 / lambda.1 #R真值
  N = 0 #循环初始化
  Y.0 = NULL
  Y.1 = NULL
  index = 0
  cat(
    'R:',
    R,
    "lam0:",
    lambda.0,
    "lam1:",
    lambda.1,
    'alpha:',
    alpha,
    "simN:",
    sim.N,
    "N0:",
    N.0,
    "N1:",
    N.1,
    '\n'
  ) #参数信息打印
  cumu2 = function(r){
  n = N.0 + N.1
  saddle.func = function(theta) {
    lam.0 * exp(theta * N.1) -  lam.1 * r * exp(-theta * N.0 * r)
  }
  theta = uniroot(saddle.func, interval = c(-5, 5), tol = 1e-15)$root
  
  phi = exp(N.0 * lam.0 * (exp(theta * N.1) - 1)) *
    exp(N.1 * lam.1 * (exp(-theta * N.0 * r) - 1)) #\varphi_Z(\theta) K_Z(\theta)
  
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
  if (is.nan(B0)) {
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
  return(SPA)}  #废案 #SPA main function
  
  SPA = function(r) {
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
    P = pnorm(w) + dnorm(v) * (1 / w - 1 / v)
    return(P)
  }

  f = function(u) {
    SPA(u) - alpha / 2
  }
  g = function(u) {
    SPA(u) - (1 - alpha / 2)
  }

    
  while (N < sim.N){
    cat(paste("模拟进度:", N/sim.N*100, "%"),'\r')
    #for(j in 1:N.0)
    #{    
    #  Y.0[j] = rpois(1, lambda.0)
    #  Y.1[j] = rpois(1, lambda.1)
    #}
    Y.0 = rpois(n = N.0, lambda = lambda.0) #sampling
    Y.1 = rpois(n = N.1, lambda = lambda.1) #sampling
    if (sum(Y.0) >= 1 & sum(Y.1) >= 1){
      lam.0 = mean(Y.0) #MLE
      lam.1 = mean(Y.1) #MLE
      N = N + 1
    }
    else{next}
    R.hat = lam.0 / lam.1 #计算R估计
#    if (lam.0 ==0 | lam.1==0){
#      next
#    }
    Rhat = c(Rhat, R.hat)
    #R.hat.var = R.hat^2*( lam.0 / (lam.0 - lam.0^2) + lam.1 / (lam.1 - lam.1^2) ) #计算R估计的方差
    #R.hat.var = R.hat^2*( lam.0 / (lam.0 - lam.0^2) + lam.1 / (lam.1 - lam.1^2) ) #计算R估计的方差
    R.hat.var = (lam.0/N.0)/lam.1^2 + (lam.0^2 * (lam.1/N.1) ) / lam.1^4
    if (length(R.hat.var) == 0 | R.hat.var <= 0 | R.hat.var == Inf){
      R.hat.var = (lam.0 / N.0)/(lam.1 / N.1)
    }
#################SPA#################
    #root1 = bisection(f, interval = c(0.0001, R.hat-0.0001), tol = 1e-5, max_iter = 1000)
    #root1 = uniroot(f, interval = c(0.00001, R.hat-0.0001), tol = 1e-15)$root
    root1 = uniroot(f, interval = c(0.01, 10), tol = 1e-15)$root

    
    if(is.null(root1)){
      index = index + 1
      cat('无解')
      fault = fault + 1
      next
    }
    #root2 = bisection(f, interval = c(R.hat+0.0001, 2*R.hat), tol = 1e-5, max_iter = 1000)
    #root2 = uniroot(f, interval = c(R.hat + 0.0001, 3*R.hat), tol = 1e-15)$root
    root2 = uniroot(g, interval = c(0.01, 10), tol = 1e-15)$root
    
    if(is.null(root2)){
      cat('无解')
      fault = fault + 1
      next
    }

    #cat('root1:',root1,'root2:',root2, "R.hat:",R.hat,'iter:',N,'\n')
    if (root1 < root2) {
      SPA.L = root1
      SPA.U = root2
    }
    else{
      SPA.U = root1
      SPA.L = root2
    }
    if (SPA.L <= R & R <= SPA.U) {
      #cat("覆盖！", "low:", SPA.L,'up:', SPA.U , 'length:', (SPA.U - SPA.L), '\n')
      SPA.K = SPA.K + 1
      #cat(SPA.K,'\n')
      SPA.Low = c(SPA.Low, SPA.L)
      SPA.Up = c(SPA.Up, SPA.U)
    }
#################Wald#################    
    Wald.L = max((R.hat - Z * sqrt(R.hat.var)), 0)
    Wald.U = R.hat + Z * sqrt(R.hat.var)
    if (Wald.L <= R & R <= Wald.U){
      Wald.K = Wald.K + 1
      Wald.Low = c(Wald.Low, Wald.L)
      Wald.Up = c(Wald.Up, Wald.U)
    }

    #cat('Wald.L:', Wald.L, '\n')
    
################Woolf#################
    Woolf.L = max(exp(log(R.hat) - Z * (1 / sum(Y.1) + 1 / (N.1 - sum(Y.1)) + 
                                        1 / sum(Y.0) + 1 / (N.0 - sum(Y.0))) ^ 0.5), 0)
    Woolf.U = exp(log(R.hat) + Z * (1 / sum(Y.1) + 1 / (N.1 - sum(Y.1)) + 1 /
                                    sum(Y.0) + 1 / (N.0 - sum(Y.0))) ^ 0.5)
    
    if (!is.nan(Woolf.L) & !is.nan(Woolf.U) & !is.infinite(Woolf.L) & !is.infinite(Woolf.U)){
      if (Woolf.L <= R & R <= Woolf.U){
        Woolf.K = Woolf.K + 1
        Woolf.Low = c(Woolf.Low, Woolf.L)
        Woolf.Up = c(Woolf.Up, Woolf.U)
      }
      #cat("low:",Woolf.L,'up:',Woolf.U,'length:',(Woolf.U-Woolf.L),'\n')
      #cat(Y.0[1],'\n')
    }
    
################LRT#################
    LRT.L = max((R.hat - sqrt(chsiq*R.hat.var)), 0)
    LRT.U = R.hat + sqrt(chsiq*R.hat.var)
    #print(cat("low:", LRT.L,'up:', LRT.U, 'length:', (LRT.U-LRT.L)))
    if (!is.nan(LRT.L) & !is.nan(LRT.U) & !is.infinite(LRT.L) & !is.infinite(LRT.U)){
      if (LRT.L <= R & R <= LRT.U){
        LRT.K = LRT.K + 1
        LRT.Low = c(LRT.Low, LRT.L)
        LRT.Up = c(LRT.Up, LRT.U)
      }
    }
################log################# 
    p.0 = lam.0
    p.1 = lam.1
    log.L = max((exp(log(p.0/p.1) - Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + 
                                              ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))), 0)
    log.U = min((exp(log(p.0/p.1) + Z*sqrt( ((N.1-sum(Y.1))/sum(Y.1))/N.1 + 
                                              ((N.0-sum(Y.0))/sum(Y.0))/N.0 ))), 3*R.hat)
    #print(cat("low:",log.L,'up:',log.U ,'length:',(log.U - log.L)))
    if (!is.nan(log.L) & !is.nan(log.U) & !is.infinite(log.L) & !is.infinite(log.U)){
      if (log.L <= R & R <= log.U){
        log.K = log.K + 1
        log.Low = c(log.Low, log.L)
        log.Up = c(log.Up, log.U)
      }
    }    
    
################Fieller#################     
    E = lam.1 ^ (2 / 3) - Z ^ 2 * (1 - lam.1) / (9 * N.1 * lam.1 ^ (1 / 3))
    D = (lam.1 * lam.0) ^ (1 / 3)
    F = lam.0 ^ (2 / 3) - Z ^ 2 * (1 - lam.0) / (9 * N.0 * lam.0 ^ (1 / 3))
    F.L = max((((D - sqrt(D ^ 2 - E * F)) / E) ^ 3), 0)
    F.U = min((((D + sqrt(D ^ 2 - E * F )) / E) ^ 3), 3 * R.hat)
    
    
    
    #a = lam.0
    #b = lam.1
    #a.var = var(Y.0)
    #b.var = var(Y.1)
    #c12 = cov(Y.0,Y.1)
    
    #q.ex = a^2/a.var
    #q.com = (b^2*a.var - 2*a*b*c12+ a^2*b.var) / (a.var*b.var-c12^2)
    #q = qt(p = 1-0.05/2, df = (length(Y.0)+length(Y.1)) )

    #F.L = 1/(a^2-q^2*a.var) * ( (a*b-q^2*c12) + sqrt( (a*b-q^2*c12)^2  - (a^2-q^2*a.var)*(b^2-q^2*b.var)  )  ) 
    
    #F.U = 1/(a^2-q^2*a.var) * ( (a*b-q^2*c12) - sqrt( (a*b-q^2*c12)^2  - (a^2-q^2*a.var)*(b^2-q^2*b.var)  )  ) 
    
    #cat("low:",F.L,'up:',F.U,'length:',(F.U-F.L))
    if (!is.nan(F.L) & !is.nan(F.U) & !is.infinite(F.L) & !is.infinite(F.U)){
      if (F.L <= R & R <= F.U){
        F.K = F.K + 1
        F.Low = c(F.Low, F.L)
        F.Up = c(F.Up, F.U)
      }
    }    
    
################Exact#################
    Exact = poisson.test(c(sum(Y.0),sum(Y.1)), T = c(N.0, N.1), r = R.hat,
                 alternative = "two.sided",
                 conf.level = 0.95)
    Exact.L = Exact$conf.int[1]
    Exact.U = Exact$conf.int[2]
    if (is.infinite(Exact.U)==TRUE){
      Exact.U = 3 * R.hat
      #cat('Exact.U:',Exact.U,'\n')
    } 
    if (!is.nan(Exact.L) & !is.nan(Exact.U) & !is.infinite(Exact.L) & !is.infinite(Exact.U)){
      if (Exact.L <= R & R <=Exact.U){
        Exact.K = Exact.K + 1
        Exact.Low = c(Exact.Low, Exact.L)
        Exact.Up = c(Exact.Up, Exact.U)
      }
    }  
  } #循环结束标志

#################结果汇总#################
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")
  cat("R真值:", R, ';', '模拟结束于:', formatted_time, '\n')
  res = data.frame(
    mean_Low = c(mean(SPA.Low),
                 mean(Wald.Low),
                 mean(Woolf.Low), 
                 mean(LRT.Low),
                 mean(log.Low),
                 mean(F.Low), 
                 mean(Exact.Low)),
    mean_Up = c(mean(SPA.Up),
                mean(Wald.Up),
                mean(Woolf.Up), 
                mean(LRT.Up), 
                mean(log.Up),
                mean(F.Up), 
                mean(Exact.Up)),
    mean_len = c(mean(SPA.Up - SPA.Low), 
                 mean(Wald.Up - Wald.Low), 
                 mean(Woolf.Up - Woolf.Low), 
                 mean(LRT.Up - LRT.Low),
                 mean(log.Up - log.Low), 
                 mean(F.Up - F.Low),
                 mean(Exact.Up - Exact.Low)),
    CR = c(SPA.K/sim.N, 
           Wald.K/sim.N,
           Woolf.K/sim.N, 
           LRT.K/sim.N, 
           log.K/sim.N, 
           F.K/sim.N,
           Exact.K/sim.N),
    mean_Rhat = mean(Rhat))
  rownames(res) <- c("SPA", "Wald", 'Woolf', 'LRT', 'log', 'Fieller', "Exact")
  return(res)
}

set.seed(NULL)
sim.N = 2000

######1.大样本,中概率,R>1######
lambda.0 = 0.8 
lambda.1 = 0.4
N.0 = 100
N.1 = 100
simulation(sim.N, N.0, N.1, lambda.0, lambda.1, alpha=0.05)


######2.小样本,小概率,R>1######
lambda.0 = 3
lambda.1 = 1.8
N.0 = 15
N.1 = 20
simulation(sim.N, N.0, N.1, lambda.0, lambda.1, 0.05)


######3.小样本,大概率,R>1######
lambda.0 = 1.2
lambda.1 = 0.2
N.0 = 1000
N.1 = 1000
simulation(sim.N, N.0, N.1, lambda.0, lambda.1, 0.05)



######4.大样本,大概率,R=1######
lambda.0 = 2
lambda.1 = 2
N.0 = 200
N.1 = 200
simulation(sim.N, N.0, N.1, lambda.0, lambda.1, 0.05)


######5.小样本,大概率,R=1######
lambda.0 = 2 
lambda.1 = 2
N.0 = 15
N.1 = 20
simulation(sim.N, N.0, N.1, lambda.0, lambda.1, 0.05)


######6.小样本,大概率,R=1######
lambda.0 = 0.75
lambda.1 = 0.75
N.0 = 200
N.1 = 200
simulation(sim.N, N.0, N.1, lambda.0, lambda.1, 0.05)


######7.大样本,小概率,R<1######
lambda.0 = 0.1
lambda.1 = 0.3
N.0 = 200
N.1 = 150
simulation(sim.N, N.0, N.1, lambda.0, lambda.1, 0.05)

######8.小样本,小概率,R<1######
lambda.0 = 1
lambda.1 = 3
N.0 = 20
N.1 = 10
simulation(sim.N, N.0, N.1, lambda.0, lambda.1, 0.05)
  

######9.小样本,小概率,R<1######
lambda.0 = 0.05
lambda.1 = 0.09
N.0 = 400
N.1 = 350
simulation(sim.N, N.0, N.1, lambda.0, lambda.1, 0.05)



