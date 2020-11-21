#Note: The calculations in this script take only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")
library("copula")

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

###################################### CoVaR calculations ##########################################
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

#VaR of BTC
VaR_rBTC = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", 0.05, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])

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

#statistics
stats = function(ts){
  vec = as.vector(ts)
  return(round(c("Min" = min(vec),
                 "Mean" = mean(vec),
                 "Median" = median(vec),
                 "Max" = max(vec),
                 "Sd" = sd(vec)),4))
}
tbl = sapply(list(VaR_rBTC, t_CoVaR_1, t_CoVaR_2, t_CoVaR_3, t_CoVaR_4, t_MCoVaR, t_VCoVaR), stats)
colnames(tbl) = c("VaR(BTC)", "CoVaR|LTC", "CoVaR|XMR", "CoVaR|XRP", "SCoVaR", "MCoVaR", "VCoVaR")
