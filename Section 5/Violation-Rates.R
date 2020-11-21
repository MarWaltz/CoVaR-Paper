#Note: The calculations in this script take approximately 1h on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")
library("rmgarch")
library("copula")
library("nloptr")

#load BTC, LTC, XMR, XRP from 'https://github.com/MarWaltz/Thesis/tree/master/Section 5/Data'
rBTC = diff(log(BTC))
rLTC = diff(log(LTC))
rXMR = diff(log(XMR))
rXRP = diff(log(XRP))
rSys = rLTC + rXMR + rXRP

set.seed(312)

################################### Return Margin Modeling ########################################
#BTC
specrBTC = ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean = 1),
                      variance.model = list(model = "gjrGARCH", garchOrder = c(3,0)),
                      distribution.model = "sstd")
fitrBTC = ugarchfit(spec = specrBTC, data = rBTC, solver = "hybrid")
v = as.vector(pit(fitrBTC))

#LTC
specrLTC = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrLTC = ugarchfit(spec = specrLTC, data = rLTC, solver = "hybrid")
u1 = as.vector(pit(fitrLTC))

#XMR
specrXMR = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,1)),
                      distribution.model = "sstd")
fitrXMR = ugarchfit(spec = specrXMR, data = rXMR, solver = "hybrid")
u2 = as.vector(pit(fitrXMR))

#XRP
specrXRP = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 1),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrXRP = ugarchfit(spec = specrXRP, data = rXRP, solver = "hybrid")
u3 = as.vector(pit(fitrXRP))

#System
specrSys = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,0)),
                      distribution.model = "sstd")
fitrSys = ugarchfit(spec = specrSys, data = rSys, solver = "hybrid")
u4 = as.vector(pit(fitrSys))

####################################### UNIVARIATE VAR ############################################
VaR_rBTC = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", 0.05, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])
VaR_rLTC = as.vector(fitted(fitrLTC)) + as.vector(sigma(fitrLTC))*qdist("sstd", 0.05, shape = coef(fitrLTC)["shape"], skew = coef(fitrLTC)["skew"])
VaR_rXMR = as.vector(fitted(fitrXMR)) + as.vector(sigma(fitrXMR))*qdist("sstd", 0.05, shape = coef(fitrXMR)["shape"], skew = coef(fitrXMR)["skew"])
VaR_rXRP = as.vector(fitted(fitrXRP)) + as.vector(sigma(fitrXRP))*qdist("sstd", 0.05, shape = coef(fitrXRP)["shape"], skew = coef(fitrXRP)["skew"])
VaR_rSys = as.vector(fitted(fitrSys)) + as.vector(sigma(fitrSys))*qdist("sstd", 0.05, shape = coef(fitrSys)["shape"], skew = coef(fitrSys)["skew"])

#BTC
VaR_BTC_exc = which(rBTC <= VaR_rBTC)
VaR_BTC_exc_rate = length(VaR_BTC_exc)/length(rBTC)

#LTC
VaR_LTC_exc = which(rLTC <= VaR_rLTC)
VaR_LTC_exc_rate = length(VaR_LTC_exc)/length(rLTC)

#XMR
VaR_XMR_exc = which(rXMR <= VaR_rXMR)
VaR_XMR_exc_rate = length(VaR_XMR_exc)/length(rXMR)

#XRP
VaR_XRP_exc = which(rXRP <= VaR_rXRP)
VaR_XRP_exc_rate = length(VaR_XRP_exc)/length(rXRP)

#System
VaR_Sys_exc = which(rSys <= VaR_rSys)
VaR_Sys_exc_rate = length(VaR_Sys_exc)/length(rSys)

########################### Copula estimates and CoVaR calculations ###############################
get_BivaCoVaR = function(copula, theta, df = NULL, cond.means, cond.sigmas, shapeY = NULL, 
                         skewY = NULL, alpha = 0.05, beta = 0.05){
  if(copula == "normal"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = normalCopula(param = theta)) - alpha*beta)
    }
  }
  if(copula == "t"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = tCopula(param = theta, df = df)) - alpha*beta)
    }
  }
  if(copula == "gumbel"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = gumbelCopula(param = theta)) - alpha*beta)
    }
  }
  if(copula == "clayton"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = claytonCopula(param = theta)) - alpha*beta)
    }
  }
  tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
  CoVaR = as.vector(cond.means) + as.vector(cond.sigmas)*qdist("sstd", tmp, shape = shapeY, skew = skewY)  
  
  return(CoVaR)
}

#################################### normal-COPULA COVARS ########################################
#BTC-LTC
TIcop1 = fitCopula(copula = normalCopula(), data = cbind(v, u1), estimate.variance = F, method = "ml")
norm_CoVaR_1 = get_BivaCoVaR(theta = TIcop1@estimate, cond.means = fitted(fitrBTC), 
                             cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"],
                             skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                             copula = "normal")

#BTC-XMR
TIcop2 = fitCopula(copula = normalCopula(), data = cbind(v, u2), estimate.variance = F, method = "ml")
norm_CoVaR_2 = get_BivaCoVaR(theta = TIcop2@estimate, cond.means = fitted(fitrBTC), 
                             cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"],
                             skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                             copula = "normal")

#BTC-XRP
TIcop3 = fitCopula(copula = normalCopula(), data = cbind(v, u3), estimate.variance = F, method = "ml")
norm_CoVaR_3 = get_BivaCoVaR(theta = TIcop3@estimate, cond.means = fitted(fitrBTC), 
                             cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"],
                             skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                             copula = "normal")

#SCoVaR
TIcop4 = fitCopula(copula = normalCopula(), data = cbind(v, u4), estimate.variance = F, method = "ml")
norm_CoVaR_4 = get_BivaCoVaR(theta = TIcop4@estimate, cond.means = fitted(fitrBTC), 
                             cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"],
                             skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                             copula = "normal")

#MCoVaR
TIfitMulti = fitCopula(copula = normalCopula(dim = 4, dispstr = "un"), estimate.variance = F,
                       data = cbind(v, u1, u2, u3), method = "ml")
copMulti = normalCopula(param = TIfitMulti@estimate, dim = 4, dispstr = "un")

MinF = function(v){
  return(pCopula(c(v, alpha, alpha, alpha), copula = copMulti) - 
           pCopula(c(1, alpha, alpha, alpha), copula = copMulti)*beta)
}
alpha = beta = 0.05; set.seed(123)
tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
norm_MCoVaR = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", tmp, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])

#VCoVaR
copVul = rotCopula(copMulti)

MinF = function(v){
  return((v - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul) + 
            pCopula(c(1-v, 1-alpha, 1-alpha, 1-alpha), copula = copVul)) -
           (beta * (1 - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul))))
}
alpha = beta = 0.05; set.seed(123)
tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
norm_VCoVaR = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", tmp, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])

###################################### t-COPULA COVARS ###########################################
#BTC-LTC
TIcop1 = fitCopula(copula = tCopula(), data = cbind(v, u1), estimate.variance = F, method = "ml")
t_CoVaR_1 = get_BivaCoVaR(theta = TIcop1@estimate[1], df = round(TIcop1@estimate[2]), 
                          cond.means = fitted(fitrBTC), cond.sigmas = sigma(fitrBTC), 
                          shapeY = coef(fitrBTC)["shape"], skewY = coef(fitrBTC)["skew"], 
                          alpha = 0.05, beta = 0.05, copula = "t")

#BTC-XMR
TIcop2 = fitCopula(copula = tCopula(), data = cbind(v, u2), estimate.variance = F, method = "ml")
t_CoVaR_2 = get_BivaCoVaR(theta = TIcop2@estimate[1], df = round(TIcop2@estimate[2]), 
                          cond.means = fitted(fitrBTC), cond.sigmas = sigma(fitrBTC), 
                          shapeY = coef(fitrBTC)["shape"], skewY = coef(fitrBTC)["skew"], 
                          alpha = 0.05, beta = 0.05, copula = "t")

#BTC-XRP
TIcop3 = fitCopula(copula = tCopula(), data = cbind(v, u3), estimate.variance = F, method = "ml")
t_CoVaR_3 = get_BivaCoVaR(theta = TIcop3@estimate[1], df = round(TIcop3@estimate[2]), 
                          cond.means = fitted(fitrBTC), cond.sigmas = sigma(fitrBTC), 
                          shapeY = coef(fitrBTC)["shape"], skewY = coef(fitrBTC)["skew"], 
                          alpha = 0.05, beta = 0.05, copula = "t")

#SCoVaR
TIcop4 = fitCopula(copula = tCopula(), data = cbind(v, u4), estimate.variance = F, method = "ml")
t_CoVaR_4 = get_BivaCoVaR(theta = TIcop4@estimate[1], df = round(TIcop4@estimate[2]), 
                          cond.means = fitted(fitrBTC), cond.sigmas = sigma(fitrBTC), 
                          shapeY = coef(fitrBTC)["shape"], skewY = coef(fitrBTC)["skew"], 
                          alpha = 0.05, beta = 0.05, copula = "t")

#MCoVaR
TIfitMulti = fitCopula(copula = tCopula(dim = 4, dispstr = "un"), estimate.variance = F,
                       data = cbind(v, u1, u2, u3), method = "ml")
copMulti = tCopula(param = TIfitMulti@estimate[-7], df = round(TIfitMulti@estimate[7]), dim = 4, dispstr = "un")

MinF = function(v){
  return(pCopula(c(v, alpha, alpha, alpha), copula = copMulti) - 
           pCopula(c(1, alpha, alpha, alpha), copula = copMulti)*beta)
}
alpha = beta = 0.05; set.seed(123)
tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
t_MCoVaR = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", tmp, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])

#VCoVaR
copVul = rotCopula(copMulti)

MinF = function(v){
  return((v - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul) + 
            pCopula(c(1-v, 1-alpha, 1-alpha, 1-alpha), copula = copVul)) -
           (beta * (1 - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul))))
}
alpha = beta = 0.05; set.seed(123)
tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
t_VCoVaR = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", tmp, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])

#################################### clayton-COPULA COVARS ########################################
#BTC-LTC
TIcop1 = fitCopula(copula = claytonCopula(), data = cbind(v, u1), estimate.variance = F, method = "ml")
Clayton_CoVaR_1 = get_BivaCoVaR(theta = TIcop1@estimate, cond.means = fitted(fitrBTC), 
                                cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"], 
                                skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                                copula = "clayton")

#BTC-XMR
TIcop2 = fitCopula(copula = claytonCopula(), data = cbind(v, u2), estimate.variance = F, method = "ml")
Clayton_CoVaR_2 = get_BivaCoVaR(theta = TIcop2@estimate, cond.means = fitted(fitrBTC), 
                                cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"], 
                                skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                                copula = "clayton")

#BTC-XRP
TIcop3 = fitCopula(copula = claytonCopula(), data = cbind(v, u3), estimate.variance = F, method = "itau") #ML estimation fails
Clayton_CoVaR_3 = get_BivaCoVaR(theta = TIcop3@estimate, cond.means = fitted(fitrBTC), 
                                cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"], 
                                skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                                copula = "clayton")

#SCoVaR
TIcop4 = fitCopula(copula = claytonCopula(), data = cbind(v, u4), estimate.variance = F, method = "ml")
Clayton_CoVaR_4 = get_BivaCoVaR(theta = TIcop4@estimate, cond.means = fitted(fitrBTC),
                                cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"], 
                                skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                                copula = "clayton")

#MCoVaR
TIfitMulti = fitCopula(copula = claytonCopula(dim = 4), estimate.variance = F,
                       data = cbind(v, u1, u2, u3), method = "ml")
copMulti = claytonCopula(param = TIfitMulti@estimate, dim = 4)

MinF = function(v){
  return(pCopula(c(v, alpha, alpha, alpha), copula = copMulti) - 
           pCopula(c(1, alpha, alpha, alpha), copula = copMulti)*beta)
}
alpha = beta = 0.05; set.seed(123)
tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
Clayton_MCoVaR = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", tmp, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])

#VCoVaR
copVul = rotCopula(copMulti)

MinF = function(v){
  return((v - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul) + 
            pCopula(c(1-v, 1-alpha, 1-alpha, 1-alpha), copula = copVul)) -
           (beta * (1 - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul))))
}
alpha = beta = 0.05; set.seed(123)
tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
Clayton_VCoVaR = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", tmp, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])

#################################### gumbel-COPULA COVARS #######################################
#BTC-LTC
TIcop1 = fitCopula(copula = gumbelCopula(), data = cbind(v, u1), estimate.variance = F, method = "ml")
Gumbel_CoVaR_1 = get_BivaCoVaR(theta = TIcop1@estimate, cond.means = fitted(fitrBTC), 
                               cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"], 
                               skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                               copula = "gumbel")

#BTC-XMR
TIcop2 = fitCopula(copula = gumbelCopula(), data = cbind(v, u2), estimate.variance = F, method = "ml")
Gumbel_CoVaR_2 = get_BivaCoVaR(theta = TIcop2@estimate, cond.means = fitted(fitrBTC), 
                               cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"], 
                               skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                               copula = "gumbel")

#BTC-XRP
TIcop3 = fitCopula(copula = gumbelCopula(), data = cbind(v, u3), estimate.variance = F, method = "ml")
Gumbel_CoVaR_3 = get_BivaCoVaR(theta = TIcop3@estimate, cond.means = fitted(fitrBTC), 
                               cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"], 
                               skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                               copula = "gumbel")

#SCoVaR
TIcop4 = fitCopula(copula = gumbelCopula(), data = cbind(v, u4), estimate.variance = F, method = "ml")
Gumbel_CoVaR_4 = get_BivaCoVaR(theta = TIcop4@estimate, cond.means = fitted(fitrBTC), 
                               cond.sigmas = sigma(fitrBTC), shapeY = coef(fitrBTC)["shape"], 
                               skewY = coef(fitrBTC)["skew"], alpha = 0.05, beta = 0.05,
                               copula = "gumbel")

#MCoVaR
TIfitMulti = fitCopula(copula = gumbelCopula(dim = 4), estimate.variance = F,
                       data = cbind(v, u1, u2, u3), method = "ml")
copMulti = gumbelCopula(param = TIfitMulti@estimate, dim = 4)

MinF = function(v){
  return(pCopula(c(v, alpha, alpha, alpha), copula = copMulti) - 
           pCopula(c(1, alpha, alpha, alpha), copula = copMulti)*beta)
}
alpha = beta = 0.05; set.seed(123)
tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
Gumbel_MCoVaR = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", tmp, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])

#VCoVaR
copVul = rotCopula(copMulti)

MinF = function(v){
  return((v - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul) + 
            pCopula(c(1-v, 1-alpha, 1-alpha, 1-alpha), copula = copVul)) -
           (beta * (1 - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul))))
}
alpha = beta = 0.05; set.seed(123)
tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
Gumbel_VCoVaR = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", tmp, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])

########################### DYNAMIC COPULA FOLLOWING PATTON(2006) #################################

Transform = function(x){
  return(tanh(x/2))
}

UpdatePatton = function(thetaInitial, a, b, c, df, data){
  theta = thetaInitial
  
  for(t in 2:nrow(data)){
    if(t <= 10){
      A = mean(qt(data[1:(t-1),1], df = df)*qt(data[1:(t-1),2], df = df))
    }else{
      A = mean(qt(data[(t-10):(t-1),1], df = df)*qt(data[(t-10):(t-1),2], df = df))
    }
    theta[t] = Transform(a + b * theta[t-1] + c * A)
  }
  return(theta)
}

LogLike = function(dataVec, theta, df){
  return(dCopula(u = dataVec, copula = tCopula(param = theta, df = df), log = TRUE))
}

ObjectiveF = function(param, data){
  a = param[1]
  b = param[2]
  c = param[3]
  df = param[4]
  
  #Initialize parameter vector
  startTheta = fitCopula(copula = tCopula(df = df, df.fixed = TRUE), 
                         data = data, estimate.variance = FALSE)@estimate[1]
  theta = rep(startTheta, nrow(data))
  
  #Update parameter vector according to Patton (2006)
  theta = UpdatePatton(thetaInitial = theta, a = a, b = b, c = c, df = df, data = data)
  
  #Calculate Log-Likelihoods
  LL = apply(cbind(data, theta), 1, function(x){LogLike(dataVec = x[1:2], theta = x[3], df = df)})
  
  #return
  return(-sum(LL))
}

set.seed(123)

### Optimization
out_BTC_LTC = nloptr(x0 = c(0.5, 0.5, 0.5, 4), eval_f = ObjectiveF, lb = c(-Inf, -Inf, -Inf, 2),
                     ub = c(Inf, Inf, Inf, 200), opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-04),
                     data = cbind(v, u1))

out_BTC_XMR = nloptr(x0 = c(0.5, 0.5, 0.5, 4), eval_f = ObjectiveF, lb = c(-Inf, -Inf, -Inf, 2),
                     ub = c(Inf, Inf, Inf, 200), opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-04),
                     data = cbind(v, u2))

out_BTC_XRP = nloptr(x0 = c(0.5, 0.5, 0.5, 4), eval_f = ObjectiveF, lb = c(-Inf, -Inf, -Inf, 2),
                     ub = c(Inf, Inf, Inf, 200), opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-04),
                     data = cbind(v, u3))

out_BTC_Sys = nloptr(x0 = c(0.5, 0.5, 0.5, 4), eval_f = ObjectiveF, lb = c(-Inf, -Inf, -Inf, 2),
                     ub = c(Inf, Inf, Inf, 200), opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-04),
                     data = cbind(v, u4))

#get correlations
getThetaPath = function(out, data){
  a = out$solution[1]
  b = out$solution[2]
  c = out$solution[3]
  df = out$solution[4]
  
  startTheta = fitCopula(copula = tCopula(df = df, df.fixed = TRUE), 
                         data = data, estimate.variance = FALSE)@estimate[1]
  theta = rep(startTheta, nrow(data))
  theta = UpdatePatton(thetaInitial = theta, a = a, b = b, c = c, df = df, data = data)
  return(theta)
}

#BTC-LTC
Patton_Theta_BTC_LTC = getThetaPath(out_BTC_LTC, data = cbind(v, u1))

#BTC-XMR
Patton_Theta_BTC_XMR = getThetaPath(out_BTC_XMR, data = cbind(v, u2))

#BTC-XRP
Patton_Theta_BTC_XRP = getThetaPath(out_BTC_XRP, data = cbind(v, u3))

#BTC-System
Patton_Theta_BTC_Sys = getThetaPath(out_BTC_Sys, data = cbind(v, u4))

#Dynamic CoVaR
get_BivaCoVaR_TV = function(rho, df, cond.means, cond.sigmas, shapeY, skewY, alpha = 0.05, beta = 0.05){
  CoVaR = c()
  
  for(t in seq_along(rho)){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = tCopula(param = rho[t], df = df)) - alpha*beta)
    }
    tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
    CoVaR[t] = as.vector(cond.means)[t] + as.vector(cond.sigmas)[t]*qdist("sstd", tmp, shape = shapeY, skew = skewY)
  }
  return(CoVaR)
}

#BTC-LTC
Patton_CoVaR_1 = get_BivaCoVaR_TV(rho = Patton_Theta_BTC_LTC, df = round(out_BTC_LTC$solution[4]), 
                                  cond.means = fitted(fitrBTC), cond.sigmas = sigma(fitrBTC), 
                                  shapeY = coef(fitrBTC)["shape"], skewY = coef(fitrBTC)["skew"], 
                                  alpha = 0.05, beta = 0.05)

#BTC-XMR
Patton_CoVaR_2 = get_BivaCoVaR_TV(rho = Patton_Theta_BTC_XMR, df = round(out_BTC_XMR$solution[4]), 
                                  cond.means = fitted(fitrBTC), cond.sigmas = sigma(fitrBTC), 
                                  shapeY = coef(fitrBTC)["shape"], skewY = coef(fitrBTC)["skew"], 
                                  alpha = 0.05, beta = 0.05)

#BTC-XRP
Patton_CoVaR_3 = get_BivaCoVaR_TV(rho = Patton_Theta_BTC_XRP, df = round(out_BTC_XRP$solution[4]), 
                                  cond.means = fitted(fitrBTC), cond.sigmas = sigma(fitrBTC), 
                                  shapeY = coef(fitrBTC)["shape"], skewY = coef(fitrBTC)["skew"], 
                                  alpha = 0.05, beta = 0.05)

#SCoVaR
Patton_CoVaR_4 = get_BivaCoVaR_TV(rho = Patton_Theta_BTC_Sys, df = round(out_BTC_Sys$solution[4]), 
                                  cond.means = fitted(fitrBTC), cond.sigmas = sigma(fitrBTC), 
                                  shapeY = coef(fitrBTC)["shape"], skewY = coef(fitrBTC)["skew"], 
                                  alpha = 0.05, beta = 0.05)

###################################### DCC t-COPULA ###############################################

set.seed(123)

#BTC-LTC
BivDCCspec1 = cgarchspec(uspec = multispec(list(specrBTC, specrLTC)), dccOrder = c(1,1),
                         distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                   transformation = "parametric"))
BivDCCfit1 = cgarchfit(spec = BivDCCspec1, data = cbind(rBTC, rLTC), solver = c("hybrid", "solnp"))

#extract necessary information and compute CoVaR
BivDCCRhos1 = sapply(BivDCCfit1@mfit$Rt, function(x){x[1,2]})
cond.meansY = as.vector(fitted(BivDCCfit1)[,1])
cond.sigmasY = as.vector(sigma(BivDCCfit1)[,1])
shapeCop = round(coef(BivDCCfit1)["[Joint]mshape"])
shapeY = coef(BivDCCfit1)["[rBTC].shape"]
skewY = coef(BivDCCfit1)["[rBTC].skew"]

DCC_CoVaR_1 = get_BivaCoVaR_TV(rho = BivDCCRhos1, df = shapeCop, cond.means = cond.meansY, 
                               cond.sigmas = cond.sigmasY, shapeY = shapeY, skewY = skewY,
                               alpha = 0.05, beta = 0.05)

#BTC-XMR
BivDCCspec2 = cgarchspec(uspec = multispec(list(specrBTC, specrXMR)), dccOrder = c(1,1),
                         distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                   transformation = "parametric"))
BivDCCfit2 = cgarchfit(spec = BivDCCspec2, data = cbind(rBTC, rXMR), solver = c("hybrid", "solnp"))

#extract necessary information and compute CoVaR
BivDCCRhos2 = sapply(BivDCCfit2@mfit$Rt, function(x){x[1,2]})
cond.meansY = as.vector(fitted(BivDCCfit2)[,1])
cond.sigmasY = as.vector(sigma(BivDCCfit2)[,1])
shapeCop = round(coef(BivDCCfit2)["[Joint]mshape"])
shapeY = coef(BivDCCfit2)["[rBTC].shape"]
skewY = coef(BivDCCfit2)["[rBTC].skew"]

DCC_CoVaR_2 = get_BivaCoVaR_TV(rho = BivDCCRhos2, df = shapeCop, cond.means = cond.meansY, 
                               cond.sigmas = cond.sigmasY, shapeY = shapeY, skewY = skewY,
                               alpha = 0.05, beta = 0.05)

#BTC-XRP
BivDCCspec3 = cgarchspec(uspec = multispec(list(specrBTC, specrXRP)), dccOrder = c(1,1),
                         distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                   transformation = "parametric"))
BivDCCfit3 = cgarchfit(spec = BivDCCspec3, data = cbind(rBTC, rXRP), solver = c("hybrid", "solnp"))

#extract necessary information and compute CoVaR
BivDCCRhos3 = sapply(BivDCCfit3@mfit$Rt, function(x){x[1,2]})
cond.meansY = as.vector(fitted(BivDCCfit3)[,1])
cond.sigmasY = as.vector(sigma(BivDCCfit3)[,1])
shapeCop = round(coef(BivDCCfit3)["[Joint]mshape"])
shapeY = coef(BivDCCfit3)["[rBTC].shape"]
skewY = coef(BivDCCfit3)["[rBTC].skew"]

DCC_CoVaR_3 = get_BivaCoVaR_TV(rho = BivDCCRhos3, df = shapeCop, cond.means = cond.meansY, 
                               cond.sigmas = cond.sigmasY, shapeY = shapeY, skewY = skewY,
                               alpha = 0.05, beta = 0.05)

#SCoVaR
BivDCCspec4 = cgarchspec(uspec = multispec(list(specrBTC, specrSys)), dccOrder = c(1,1),
                         distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                   transformation = "parametric"))
BivDCCfit4 = cgarchfit(spec = BivDCCspec4, data = cbind(rBTC, rSys), solver = c("hybrid", "solnp"))

#extract necessary information and compute CoVaR
BivDCCRhos4 = sapply(BivDCCfit4@mfit$Rt, function(x){x[1,2]})
cond.meansY = as.vector(fitted(BivDCCfit4)[,1])
cond.sigmasY = as.vector(sigma(BivDCCfit4)[,1])
shapeCop = round(coef(BivDCCfit4)["[Joint]mshape"])
shapeY = coef(BivDCCfit4)["[rBTC].shape"]
skewY = coef(BivDCCfit4)["[rBTC].skew"]

DCC_CoVaR_4 = get_BivaCoVaR_TV(rho = BivDCCRhos4, df = shapeCop, cond.means = cond.meansY, 
                               cond.sigmas = cond.sigmasY, shapeY = shapeY, skewY = skewY,
                               alpha = 0.05, beta = 0.05)

#MCoVaR
DCCspec1 = cgarchspec(uspec = multispec(list(specrBTC, specrLTC, specrXMR, specrXRP)), dccOrder = c(1,1),
                      distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                transformation = "parametric"))
DCCfit1 = cgarchfit(spec = DCCspec1, data = cbind(rBTC, rLTC, rXMR, rXRP), solver = c("hybrid", "solnp"))

getDCCMultiCoVaR = function(DCCfit1, alpha = 0.05, beta = 0.05){
  n = length(DCCfit1@mfit$Rt)
  df = round(coef(DCCfit1)["[Joint]mshape"])
  shapeY = coef(DCCfit1)["[rBTC].shape"]
  skewY = coef(DCCfit1)["[rBTC].skew"]
  cond.meansY = as.vector(fitted(DCCfit1)[,1])
  cond.sigmasY = as.vector(sigma(DCCfit1)[,1])
  
  CoVaR = c()
  
  for(t in seq_len(n)){
    set.seed(t)
    cop = tCopula(param = P2p(DCCfit1@mfit$Rt[[t]]), dim = 4, df = df, df.fixed = TRUE, dispstr = "un")
    
    MinF = function(v){
      return(pCopula(c(v, alpha, alpha, alpha), copula = cop) - 
               pCopula(c(1, alpha, alpha, alpha), copula = cop)*beta)
    }
    tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
    CoVaR[t] = cond.meansY[t] + cond.sigmasY[t]*qdist("sstd", tmp, shape = shapeY, skew = skewY)
  }
  return(CoVaR)
}

DCC_MCoVaR = getDCCMultiCoVaR(DCCfit1, alpha = 0.05, beta = 0.05)

#VCoVaR
getDCCVCoVaR = function(DCCfit1, alpha = 0.05, beta = 0.05){
  n = length(DCCfit1@mfit$Rt)
  df = round(coef(DCCfit1)["[Joint]mshape"])
  shapeY = coef(DCCfit1)["[rBTC].shape"]
  skewY = coef(DCCfit1)["[rBTC].skew"]
  cond.meansY = as.vector(fitted(DCCfit1)[,1])
  cond.sigmasY = as.vector(sigma(DCCfit1)[,1])
  
  CoVaR = c()
  
  for(t in seq_len(n)){
    set.seed(t)
    cop = rotCopula(tCopula(param = P2p(DCCfit1@mfit$Rt[[t]]), dim = 4, df = df, df.fixed = TRUE, dispstr = "un"))
    
    MinF = function(v){
      return((v - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = cop) + 
                pCopula(c(1-v, 1-alpha, 1-alpha, 1-alpha), copula = cop)) -
               (beta * (1 - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = cop))))
    }
    tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
    CoVaR[t] = cond.meansY[t] + cond.sigmasY[t]*qdist("sstd", tmp, shape = shapeY, skew = skewY)
    cat(t)
  }
  return(CoVaR)
}

DCC_VCoVaR = getDCCVCoVaR(DCCfit1, alpha = 0.05, beta = 0.05)

######################################### EVALUATION ############################################
evalBivariate = function(retX, retY, VaRX, CoVaR){
  exc = which(retX <= VaRX)
  
  red_CoVaR = rep(NA, length(CoVaR))
  red_CoVaR[exc] = CoVaR[exc]
  CoVaR_exc = which(retY <= red_CoVaR)
  CoVaR_rate = length(CoVaR_exc)/length(exc)
  
  return(list("rate" = round(CoVaR_rate, 4), "abs" = length(CoVaR_exc), "CoVaR_exc" = CoVaR_exc))
}

evalMulti3X = function(retX1, retX2, retX3, retY, VaRX1, VaRX2, VaRX3, MCoVaR){
  exc1 = which(retX1 <= VaRX1)
  exc2 = which(retX2 <= VaRX2)
  exc3 = which(retX3 <= VaRX3)
  exc = sort(intersect(exc1, intersect(exc2, exc3)))
  
  red_MCoVaR = rep(NA, length(MCoVaR))
  red_MCoVaR[exc] = MCoVaR[exc]
  MCoVaR_exc = which(retY <= red_MCoVaR)
  MCoVaR_rate = length(MCoVaR_exc)/length(exc)
  
  return(list("rate" = round(MCoVaR_rate, 4), "abs" = length(MCoVaR_exc), "MCoVaR_exc" = MCoVaR_exc))
}

evalVulnerability3X = function(retX1, retX2, retX3, retY, VaRX1, VaRX2, VaRX3, VCoVaR){
  exc1 = which(retX1 <= VaRX1)
  exc2 = which(retX2 <= VaRX2)
  exc3 = which(retX3 <= VaRX3)
  exc = sort(unique(c(exc1, exc2, exc3)))
  
  red_VCoVaR = rep(NA, length(VCoVaR))
  red_VCoVaR[exc] = VCoVaR[exc]
  VCoVaR_exc = which(retY <= red_VCoVaR)
  VCoVaR_rate = length(VCoVaR_exc)/length(exc)
  
  return(list("rate" = round(VCoVaR_rate, 4), "abs" = length(VCoVaR_exc), "VCoVaR_exc" = VCoVaR_exc))
}

##################################### normal-COPULA COVAR #######################################

###Bivariate CoVaR: BTC-LTC
norm_CoVaR_LTC_exc = evalBivariate(retX = rLTC, retY = rBTC, VaRX = VaR_rLTC,
                                   CoVaR = norm_CoVaR_1)
###Bivariate CoVaR: BTC-XMR
norm_CoVaR_XMR_exc = evalBivariate(retX = rXMR, retY = rBTC, VaRX = VaR_rXMR,
                                   CoVaR = norm_CoVaR_2)
###Bivariate CoVaR: BTC-XRP
norm_CoVaR_XRP_exc = evalBivariate(retX = rXRP, retY = rBTC, VaRX = VaR_rXRP,
                                   CoVaR = norm_CoVaR_3)
###System-CoVaR
norm_CoVaR_Sys_exc = evalBivariate(retX = rSys, retY = rBTC, VaRX = VaR_rSys,
                                   CoVaR = norm_CoVaR_4)
###Multi-CoVaR
norm_CoVaR_Multi_exc = evalMulti3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC,
                                   VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                   MCoVaR = norm_MCoVaR)
###Vulnerability-CoVaR
norm_CoVaR_V_exc = evalVulnerability3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC, 
                                       VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                       VCoVaR = norm_VCoVaR)

######################################## t-COPULA COVAR ########################################

###Bivariate CoVaR: BTC-LTC
t_CoVaR_LTC_exc = evalBivariate(retX = rLTC, retY = rBTC, VaRX = VaR_rLTC,
                                CoVaR = t_CoVaR_1)
###Bivariate CoVaR: BTC-XMR
t_CoVaR_XMR_exc = evalBivariate(retX = rXMR, retY = rBTC, VaRX = VaR_rXMR,
                                CoVaR = t_CoVaR_2)
###Bivariate CoVaR: BTC-XRP
t_CoVaR_XRP_exc = evalBivariate(retX = rXRP, retY = rBTC, VaRX = VaR_rXRP,
                                CoVaR = t_CoVaR_3)
###System-CoVaR
t_CoVaR_Sys_exc = evalBivariate(retX = rSys, retY = rBTC, VaRX = VaR_rSys,
                                CoVaR = t_CoVaR_4)
###Multi-CoVaR
t_CoVaR_Multi_exc = evalMulti3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC, 
                                VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                MCoVaR = t_MCoVaR)
###Vulnerability-CoVaR
t_CoVaR_V_exc = evalVulnerability3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC, 
                                    VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                    VCoVaR = t_VCoVaR)

###################################### PATTON t-COPULA COVAR ######################################

###Bivariate CoVaR: BTC-LTC
Patton_CoVaR_LTC_exc = evalBivariate(retX = rLTC, retY = rBTC, VaRX = VaR_rLTC,
                                     CoVaR = Patton_CoVaR_1)
###Bivariate CoVaR: BTC-XMR
Patton_CoVaR_XMR_exc = evalBivariate(retX = rXMR, retY = rBTC, VaRX = VaR_rXMR,
                                     CoVaR = Patton_CoVaR_2)
###Bivariate CoVaR: BTC-XRP
Patton_CoVaR_XRP_exc = evalBivariate(retX = rXRP, retY = rBTC, VaRX = VaR_rXRP,
                                     CoVaR = Patton_CoVaR_3)
###System-CoVaR
Patton_CoVaR_Sys_exc = evalBivariate(retX = rSys, retY = rBTC, VaRX = VaR_rSys,
                                     CoVaR = Patton_CoVaR_4)

###################################### DCC t-COPULA COVAR #########################################

###Bivariate CoVaR: BTC-LTC
DCC_CoVaR_LTC_exc = evalBivariate(retX = rLTC, retY = rBTC, VaRX = VaR_rLTC,
                                  CoVaR = DCC_CoVaR_1)
###Bivariate CoVaR: BTC-XMR
DCC_CoVaR_XMR_exc = evalBivariate(retX = rXMR, retY = rBTC, VaRX = VaR_rXMR,
                                  CoVaR = DCC_CoVaR_2)
###Bivariate CoVaR: BTC-XRP
DCC_CoVaR_XRP_exc = evalBivariate(retX = rXRP, retY = rBTC, VaRX = VaR_rXRP,
                                  CoVaR = DCC_CoVaR_3)
###System-CoVaR
DCC_CoVaR_Sys_exc = evalBivariate(retX = rSys, retY = rBTC, VaRX = VaR_rSys,
                                  CoVaR = DCC_CoVaR_4)
###Multi-CoVaR
DCC_CoVaR_Multi_exc = evalMulti3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC,
                                  VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                  MCoVaR = DCC_MCoVaR)
###Vulnerability-CoVaR
DCC_CoVaR_V_exc = evalVulnerability3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC, 
                                      VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                      VCoVaR = DCC_VCoVaR)

##################################### CLAYTON-COPULA COVAR #######################################

###Bivariate CoVaR: BTC-LTC
Clayton_CoVaR_LTC_exc = evalBivariate(retX = rLTC, retY = rBTC, VaRX = VaR_rLTC,
                                      CoVaR = Clayton_CoVaR_1)
###Bivariate CoVaR: BTC-XMR
Clayton_CoVaR_XMR_exc = evalBivariate(retX = rXMR, retY = rBTC, VaRX = VaR_rXMR,
                                      CoVaR = Clayton_CoVaR_2)
###Bivariate CoVaR: BTC-XRP
Clayton_CoVaR_XRP_exc = evalBivariate(retX = rXRP, retY = rBTC, VaRX = VaR_rXRP,
                                      CoVaR = Clayton_CoVaR_3)
###System-CoVaR
Clayton_CoVaR_Sys_exc = evalBivariate(retX = rSys, retY = rBTC, VaRX = VaR_rSys,
                                      CoVaR = Clayton_CoVaR_4)
###Multi-CoVaR
Clayton_CoVaR_Multi_exc = evalMulti3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC,
                                      VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                      MCoVaR = Clayton_MCoVaR)
###Vulnerability-CoVaR
Clayton_CoVaR_V_exc = evalVulnerability3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC,
                                          VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                          VCoVaR = Clayton_VCoVaR)

##################################### GUMBEL-COPULA COVAR ######################################

###Bivariate CoVaR: BTC-LTC
Gumbel_CoVaR_LTC_exc = evalBivariate(retX = rLTC, retY = rBTC, VaRX = VaR_rLTC,
                                     CoVaR = Gumbel_CoVaR_1)
###Bivariate CoVaR: BTC-XMR
Gumbel_CoVaR_XMR_exc = evalBivariate(retX = rXMR, retY = rBTC, VaRX = VaR_rXMR,
                                     CoVaR = Gumbel_CoVaR_2)
###Bivariate CoVaR: BTC-XRP
Gumbel_CoVaR_XRP_exc = evalBivariate(retX = rXRP, retY = rBTC, VaRX = VaR_rXRP,
                                     CoVaR = Gumbel_CoVaR_3)
###System-CoVaR
Gumbel_CoVaR_Sys_exc = evalBivariate(retX = rSys, retY = rBTC, VaRX = VaR_rSys,
                                     CoVaR = Gumbel_CoVaR_4)
###Multi-CoVaR
Gumbel_CoVaR_Multi_exc = evalMulti3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC,
                                     VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                     MCoVaR = Gumbel_MCoVaR)
###Vulnerability-CoVaR
Gumbel_CoVaR_V_exc = evalVulnerability3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC,
                                         VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                         VCoVaR = Gumbel_VCoVaR)

############################# COMPARISON OF RELATIVE EXCEEDENCES ##################################

norm_CoVaR_rates = c(norm_CoVaR_LTC_exc[[1]], norm_CoVaR_XMR_exc[[1]],
                     norm_CoVaR_XRP_exc[[1]], norm_CoVaR_Sys_exc[[1]],
                     norm_CoVaR_Multi_exc[[1]], norm_CoVaR_V_exc[[1]])
t_CoVaR_rates = c(t_CoVaR_LTC_exc[[1]], t_CoVaR_XMR_exc[[1]],
                  t_CoVaR_XRP_exc[[1]], t_CoVaR_Sys_exc[[1]],
                  t_CoVaR_Multi_exc[[1]], t_CoVaR_V_exc[[1]])
Patton_t_CoVaR_rates = c(Patton_CoVaR_LTC_exc[[1]], Patton_CoVaR_XMR_exc[[1]],
                         Patton_CoVaR_XRP_exc[[1]], Patton_CoVaR_Sys_exc[[1]], NA, NA)
DCC_t_CoVaR_rates = c(DCC_CoVaR_LTC_exc[[1]], DCC_CoVaR_XMR_exc[[1]], DCC_CoVaR_XRP_exc[[1]],
                      DCC_CoVaR_Sys_exc[[1]], DCC_CoVaR_Multi_exc[[1]], DCC_CoVaR_V_exc[[1]])
Clayton_CoVaR_rates = c(Clayton_CoVaR_LTC_exc[[1]], Clayton_CoVaR_XMR_exc[[1]],
                        Clayton_CoVaR_XRP_exc[[1]], Clayton_CoVaR_Sys_exc[[1]],
                        Clayton_CoVaR_Multi_exc[[1]], Clayton_CoVaR_V_exc[[1]])
Gumbel_CoVaR_rates = c(Gumbel_CoVaR_LTC_exc[[1]], Gumbel_CoVaR_XMR_exc[[1]],
                       Gumbel_CoVaR_XRP_exc[[1]], Gumbel_CoVaR_Sys_exc[[1]],
                       Gumbel_CoVaR_Multi_exc[[1]], Gumbel_CoVaR_V_exc[[1]])

Exc_rates = cbind(norm_CoVaR_rates, t_CoVaR_rates, Patton_t_CoVaR_rates, 
                  DCC_t_CoVaR_rates, Clayton_CoVaR_rates, Gumbel_CoVaR_rates)
rownames(Exc_rates) = c("Bivariate: LTC", "Bivariate: XMR", "Bivariate: XRP", "System", "Multi", "Vulnerability")

############################## COMPARISON OF ABSOLUTE EXCEEDENCES #################################

CoVar_theo_Exc = c("Bivariate: LTC" = length(VaR_LTC_exc), "Bivariate: XMR" = length(VaR_XMR_exc), 
                   "Bivariate: XRP" = length(VaR_XRP_exc), "System" = length(VaR_Sys_exc), 
                   "Multi" = length(intersect(VaR_LTC_exc, intersect(VaR_XMR_exc, VaR_XRP_exc))),
                   "Vulnerability" = length(unique(c(VaR_LTC_exc, VaR_XMR_exc, VaR_XRP_exc))))*0.05
norm_CoVaR_abs = c(norm_CoVaR_LTC_exc[[2]], norm_CoVaR_XMR_exc[[2]],
                   norm_CoVaR_XRP_exc[[2]], norm_CoVaR_Sys_exc[[2]],
                   norm_CoVaR_Multi_exc[[2]], norm_CoVaR_V_exc[[2]])
t_CoVaR_abs = c(t_CoVaR_LTC_exc[[2]], t_CoVaR_XMR_exc[[2]],
                t_CoVaR_XRP_exc[[2]], t_CoVaR_Sys_exc[[2]],
                t_CoVaR_Multi_exc[[2]], t_CoVaR_V_exc[[2]])
Patton_t_CoVaR_abs = c(Patton_CoVaR_LTC_exc[[2]], Patton_CoVaR_XMR_exc[[2]],
                       Patton_CoVaR_XRP_exc[[2]], Patton_CoVaR_Sys_exc[[2]], NA, NA)
DCC_t_CoVaR_abs = c(DCC_CoVaR_LTC_exc[[2]], DCC_CoVaR_XMR_exc[[2]], DCC_CoVaR_XRP_exc[[2]],
                    DCC_CoVaR_Sys_exc[[2]], DCC_CoVaR_Multi_exc[[2]], DCC_CoVaR_V_exc[[2]])
Clayton_CoVaR_abs = c(Clayton_CoVaR_LTC_exc[[2]], Clayton_CoVaR_XMR_exc[[2]],
                      Clayton_CoVaR_XRP_exc[[2]], Clayton_CoVaR_Sys_exc[[2]],
                      Clayton_CoVaR_Multi_exc[[2]], Clayton_CoVaR_V_exc[[2]])
Gumbel_CoVaR_abs = c(Gumbel_CoVaR_LTC_exc[[2]], Gumbel_CoVaR_XMR_exc[[2]],
                     Gumbel_CoVaR_XRP_exc[[2]], Gumbel_CoVaR_Sys_exc[[2]],
                     Gumbel_CoVaR_Multi_exc[[2]], Gumbel_CoVaR_V_exc[[2]])

Exc_abs = cbind(norm_CoVaR_abs, t_CoVaR_abs, Patton_t_CoVaR_abs, 
                DCC_t_CoVaR_abs, Clayton_CoVaR_abs, Gumbel_CoVaR_abs)
rownames(Exc_abs) = c("Bivariate: LTC", "Bivariate: XMR", "Bivariate: XRP", "System", "Multi", "Vulnerability")

######################################## FINAL RESULTS ############################################
round(Exc_rates,4)
Exc_abs
