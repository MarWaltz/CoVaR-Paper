#Note: The calculations in this script take approximately 30h on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")
library("copula")

#load BTC, LTC, XMR, XRP from 'https://github.com/MarWaltz/Thesis/tree/master/Section 5/Data'
rBTC = diff(log(BTC))
rLTC = diff(log(LTC))
rXMR = diff(log(XMR))
rXRP = diff(log(XRP))
rSys = rLTC + rXMR + rXRP

#functions
get_BivaCoVaR = function(copula, theta, df = NULL, foreY, alpha = 0.05, beta = 0.05){
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
  if(copula == "clayton"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = claytonCopula(param = theta)) - alpha*beta)
    }
  }
  if(copula == "gumbel"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = gumbelCopula(param = theta)) - alpha*beta)
    }
  }
  tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
  CoVaR = as.vector(quantile(foreY, tmp))
  
  return(CoVaR)
}

get_MCoVaR = function(copula, theta, df = NULL, foreY, alpha = 0.05, beta = 0.05){
  if(copula == "normal"){
    copMulti = normalCopula(param = theta, dim = 4, dispstr = "un")
  }
  if(copula == "t"){
    copMulti = tCopula(param = theta, df = df, dim = 4, dispstr = "un")
  }
  if(copula == "clayton"){
    copMulti = claytonCopula(param = theta, dim = 4)
  }
  if(copula == "gumbel"){
    copMulti = gumbelCopula(param = theta, dim = 4)
  }
  MinF = function(v){
    return(pCopula(c(v, alpha, alpha, alpha), copula = copMulti) - 
             pCopula(c(1, alpha, alpha, alpha), copula = copMulti)*beta)
  }
  tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
  MCoVaR = as.vector(quantile(foreY, tmp))
  return(MCoVaR)
}

get_VCoVaR = function(copula, theta, df = NULL, foreY, alpha = 0.05, beta = 0.05){
  if(copula == "normal"){
    copVul = rotCopula(normalCopula(param = theta, dim = 4, dispstr = "un"))
  }
  if(copula == "t"){
    copVul = rotCopula(tCopula(param = theta, df = df, dim = 4, dispstr = "un"))
  }
  if(copula == "clayton"){
    copVul = rotCopula(claytonCopula(param = theta, dim = 4))
  }
  if(copula == "gumbel"){
    copVul = rotCopula(gumbelCopula(param = theta, dim = 4))
  }
  MinF = function(v){
    return((v - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul) + 
              pCopula(c(1-v, 1-alpha, 1-alpha, 1-alpha), copula = copVul)) -
             (beta * (1 - pCopula(c(1, 1-alpha, 1-alpha, 1-alpha), copula = copVul))))
  }
  tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
  VCoVaR = as.vector(quantile(foreY, tmp))
  return(VCoVaR)
}

########################### OUT OF SAMPLE FORECASTING USING ROLLING WINDOW ########################

getVaR_OOS = function(CurrSample, n.ahead, alpha){
  #Margin model estimation
  spec = ugarchspec(variance.model = list(garchOrder = c(1,1), model = "gjrGARCH"),
                    mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
                    distribution.model = "sstd")
  fit = ugarchfit(spec, data = CurrSample, solver = "hybrid")
  
  #forecast margin
  fore = ugarchforecast(fitORspec = fit, n.ahead = n.ahead)
  
  #calculate VaR
  FC = as.vector(quantile(fore, alpha))
  return(FC)
}


getCoVaR_OOS = function(CurrSample, type, copula, n.ahead, alpha, beta){
  if(!is.element(type, c("Bivariate", "Multi", "Vulnerability"))){
    stop("Please select valide type.")
  }
  CurrSample = as.matrix(CurrSample)
  
  #Margin model estimation
  spec = ugarchspec(variance.model = list(garchOrder = c(1,1), model = "gjrGARCH"),
                    mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
                    distribution.model = "sstd")
  fitY = ugarchfit(spec, data = CurrSample[,1], solver = "hybrid")
  v = as.vector(pit(fitY))
  
  if(type == "Bivariate"){
    fitX1 = ugarchfit(spec, data = CurrSample[,2], solver = "hybrid")
    u1 = as.vector(pit(fitX1))
  }
  if((type == "Multi") | (type == "Vulnerability")){
    fitX1 = ugarchfit(spec, data = CurrSample[,2], solver = "hybrid")
    fitX2 = ugarchfit(spec, data = CurrSample[,3], solver = "hybrid")
    fitX3 = ugarchfit(spec, data = CurrSample[,4], solver = "hybrid")
    
    u1 = as.vector(pit(fitX1))
    u2 = as.vector(pit(fitX2))
    u3 = as.vector(pit(fitX3))
  }
  
  #copula estimation
  if(type == "Bivariate"){
    if(copula == "normal"){
      estCop = try(fitCopula(copula = normalCopula(), data = cbind(v, u1), 
                             estimate.variance = F, method = "ml")@estimate, silent = T)
      if(inherits(estCop, "try-error")){
        estCop = fitCopula(copula = normalCopula(), data = cbind(v, u1), 
                           estimate.variance = F, method = "itau")@estimate
      }
    }
    if(copula == "t"){
      estCop = try(fitCopula(copula = tCopula(), data = cbind(v, u1), 
                             estimate.variance = F, method = "ml")@estimate, silent = T)
      if(inherits(estCop, "try-error")){
        estCop = fitCopula(copula = tCopula(), data = cbind(v, u1), 
                           estimate.variance = F, method = "itau")@estimate
      }
    }
    if(copula == "clayton"){
      estCop = try(fitCopula(copula = claytonCopula(), data = cbind(v, u1), 
                             estimate.variance = F, method = "ml")@estimate, silent = T)
      if(inherits(estCop, "try-error")){
        estCop = fitCopula(copula = claytonCopula(), data = cbind(v, u1), 
                           estimate.variance = F, method = "itau")@estimate
      }
    }
    if(copula == "gumbel"){
      estCop = try(fitCopula(copula = gumbelCopula(), data = cbind(v, u1), 
                             estimate.variance = F, method = "ml")@estimate, silent = T)
      if(inherits(estCop, "try-error")){
        estCop = fitCopula(copula = gumbelCopula(), data = cbind(v, u1), 
                           estimate.variance = F, method = "itau")@estimate
      }
    }
  }
  
  if((type == "Multi") | (type == "Vulnerability")){
    if(copula == "normal"){
      estCop = try(fitCopula(copula = normalCopula(dim = 4, dispstr = "un"), data = cbind(v, u1, u2, u3), 
                             estimate.variance = F, method = "ml")@estimate, silent = T)
      if(inherits(estCop, "try-error")){
        estCop = fitCopula(copula = normalCopula(dim = 4, dispstr = "un"), data = cbind(v, u1, u2, u3), 
                           estimate.variance = F, method = "itau")@estimate
      }
    }
    if(copula == "t"){
      estCop = try(fitCopula(copula = tCopula(dim = 4, dispstr = "un"), data = cbind(v, u1, u2, u3), 
                             estimate.variance = F, method = "ml")@estimate, silent = T)
      if(inherits(estCop, "try-error")){
        estCop = fitCopula(copula = tCopula(dim = 4, dispstr = "un"), data = cbind(v, u1, u2, u3), 
                           estimate.variance = F, method = "itau")@estimate
      }
    }
    if(copula == "clayton"){
      estCop = try(fitCopula(copula = claytonCopula(dim = 4), data = cbind(v, u1, u2, u3), 
                             estimate.variance = F, method = "ml")@estimate, silent = T)
      if(inherits(estCop, "try-error")){
        estCop = fitCopula(copula = claytonCopula(dim = 4), data = cbind(v, u1, u2, u3), 
                           estimate.variance = F, method = "itau")@estimate
      }
    }
    if(copula == "gumbel"){
      estCop = try(fitCopula(copula = gumbelCopula(dim = 4), data = cbind(v, u1, u2, u3), 
                             estimate.variance = F, method = "ml")@estimate, silent = T)
      if(inherits(estCop, "try-error")){
        estCop = fitCopula(copula = gumbelCopula(dim = 4), data = cbind(v, u1, u2, u3), 
                           estimate.variance = F, method = "itau")@estimate
      }
    }
  }
  
  #forecast margin of Y
  foreY = ugarchforecast(fitORspec = fitY, n.ahead = n.ahead)
  
  #forecast CoVaR measure
  if(type == "Bivariate"){
    FC = get_BivaCoVaR(copula = copula, theta = estCop[1], df = round(estCop[2]),
                       foreY = foreY, alpha = alpha, beta = beta)
  }
  if(type == "Multi"){
    if(copula == "t"){
      FC = get_MCoVaR(copula = copula, theta = estCop[-7], df = round(estCop[7]),
                      foreY = foreY, alpha = alpha, beta = beta)
    }
    if(is.element(copula, c("normal", "clayton", "gumbel"))){
      FC = get_MCoVaR(copula = copula, theta = estCop, foreY = foreY, 
                      alpha = alpha, beta = beta)
    }
  }
  if(type == "Vulnerability"){
    if(copula == "t"){
      FC = get_VCoVaR(copula = copula, theta = estCop[-7], df = round(estCop[7]),
                      foreY = foreY, alpha = alpha, beta = beta)
    }
    if(is.element(copula, c("normal", "clayton", "gumbel"))){
      FC = get_VCoVaR(copula = copula, theta = estCop, foreY = foreY, 
                      alpha = alpha, beta = beta)
    }
  }
  return(FC)
}

get_FC = function(data, type, copula, window.length = 500, n.ahead = 1, alpha = 0.05, beta = 0.05){
  #initializing
  res = c()
  data = as.matrix(data)
  nobs = nrow(data)
  stop = FALSE; set.seed(123); i = 1
  
  #run through data
  repeat{
    SampInd = (1+n.ahead*(i-1)):(window.length + n.ahead*(i-1))
    if(any((SampInd + n.ahead) > nobs)){
      if(n.ahead == 1){
        break
      }else{
        n.ahead = n.ahead - (max(SampInd + n.ahead) - nobs)
        stop = TRUE
      }
    }
    if(type == "VaR"){
      tmp = try(getVaR_OOS(CurrSample = data[SampInd,], n.ahead = n.ahead, alpha = alpha))
    }else{
      tmp = try(getCoVaR_OOS(CurrSample = data[SampInd,], n.ahead = n.ahead,
                             copula = copula, type = type, alpha = alpha, beta = beta))
    }
    if(inherits(tmp, "try-error")){
      tmp = NA
    }
    res = c(res, tmp)
    cat(i)
    if(stop == TRUE){break}else{i = i + 1}
  }
  
  #return forecasts
  return(res)
}

#out-of-sample VaR
VaRBTC = get_FC(data = rBTC, type = "VaR", window.length = 500, n.ahead = 1,
                copula = NULL, alpha = 0.05, beta = NULL)
VaRLTC = get_FC(data = rLTC, type = "VaR", window.length = 500, n.ahead = 1,
                copula = NULL, alpha = 0.05, beta = NULL)
VaRXMR = get_FC(data = rXMR, type = "VaR", window.length = 500, n.ahead = 1,
                copula = NULL, alpha = 0.05, beta = NULL)
VaRXRP = get_FC(data = rXRP, type = "VaR", window.length = 500, n.ahead = 1,
                copula = NULL, alpha = 0.05, beta = NULL)
VaRSys = get_FC(data = rSys, type = "VaR", window.length = 500, n.ahead = 1,
                copula = NULL, alpha = 0.05, beta = NULL)

#out-of-sample bivariate CoVaR and SCoVaR
norm_CoVaR_LTC = get_FC(data = cbind(rBTC, rLTC), type = "Bivariate", window.length = 500, n.ahead = 1,
                        copula = "normal", alpha = 0.05, beta = 0.05)
norm_CoVaR_XMR = get_FC(data = cbind(rBTC, rXMR), type = "Bivariate", window.length = 500, n.ahead = 1,
                        copula = "normal", alpha = 0.05, beta = 0.05)
norm_CoVaR_XRP = get_FC(data = cbind(rBTC, rXRP), type = "Bivariate", window.length = 500, n.ahead = 1,
                        copula = "normal", alpha = 0.05, beta = 0.05)
norm_SCoVaR = get_FC(data = cbind(rBTC, rSys), type = "Bivariate", window.length = 500, n.ahead = 1,
                     copula = "normal", alpha = 0.05, beta = 0.05)

t_CoVaR_LTC = get_FC(data = cbind(rBTC, rLTC), type = "Bivariate", window.length = 500, n.ahead = 1,
                     copula = "t", alpha = 0.05, beta = 0.05)
t_CoVaR_XMR = get_FC(data = cbind(rBTC, rXMR), type = "Bivariate", window.length = 500, n.ahead = 1,
                     copula = "t", alpha = 0.05, beta = 0.05)
t_CoVaR_XRP = get_FC(data = cbind(rBTC, rXRP), type = "Bivariate", window.length = 500, n.ahead = 1,
                     copula = "t", alpha = 0.05, beta = 0.05)
t_SCoVaR = get_FC(data = cbind(rBTC, rSys), type = "Bivariate", window.length = 500, n.ahead = 1,
                  copula = "t", alpha = 0.05, beta = 0.05)

Clayton_CoVaR_LTC = get_FC(data = cbind(rBTC, rLTC), type = "Bivariate", window.length = 500, n.ahead = 1,
                           copula = "clayton", alpha = 0.05, beta = 0.05)
Clayton_CoVaR_XMR = get_FC(data = cbind(rBTC, rXMR), type = "Bivariate", window.length = 500, n.ahead = 1,
                           copula = "clayton", alpha = 0.05, beta = 0.05)
Clayton_CoVaR_XRP = get_FC(data = cbind(rBTC, rXRP), type = "Bivariate", window.length = 500, n.ahead = 1,
                           copula = "clayton", alpha = 0.05, beta = 0.05)
Clayton_SCoVaR = get_FC(data = cbind(rBTC, rSys), type = "Bivariate", window.length = 500, n.ahead = 1,
                        copula = "clayton", alpha = 0.05, beta = 0.05)

Gumbel_CoVaR_LTC = get_FC(data = cbind(rBTC, rLTC), type = "Bivariate", window.length = 500, n.ahead = 1,
                          copula = "gumbel", alpha = 0.05, beta = 0.05)
Gumbel_CoVaR_XMR = get_FC(data = cbind(rBTC, rXMR), type = "Bivariate", window.length = 500, n.ahead = 1,
                          copula = "gumbel", alpha = 0.05, beta = 0.05)
Gumbel_CoVaR_XRP = get_FC(data = cbind(rBTC, rXRP), type = "Bivariate", window.length = 500, n.ahead = 1,
                          copula = "gumbel", alpha = 0.05, beta = 0.05)
Gumbel_SCoVaR = get_FC(data = cbind(rBTC, rSys), type = "Bivariate", window.length = 500, n.ahead = 1,
                       copula = "gumbel", alpha = 0.05, beta = 0.05)

#out-of-sample MCoVaR
norm_MCoVaR = get_FC(data = cbind(rBTC, rLTC, rXMR, rXRP), type = "Multi", window.length = 500, 
                     n.ahead = 1, copula = "normal", alpha = 0.05, beta = 0.05)
t_MCoVaR = get_FC(data = cbind(rBTC, rLTC, rXMR, rXRP), type = "Multi", window.length = 500, 
                  n.ahead = 1, copula = "t", alpha = 0.05, beta = 0.05)
Clayton_MCoVaR = get_FC(data = cbind(rBTC, rLTC, rXMR, rXRP), type = "Multi", window.length = 500, 
                        n.ahead = 1, copula = "clayton", alpha = 0.05, beta = 0.05)
Gumbel_MCoVaR = get_FC(data = cbind(rBTC, rLTC, rXMR, rXRP), type = "Multi", window.length = 500, 
                        n.ahead = 1, copula = "gumbel", alpha = 0.05, beta = 0.05)

#out-of-sample VCoVaR
norm_VCoVaR = get_FC(data = cbind(rBTC, rLTC, rXMR, rXRP), type = "Vulnerability", window.length = 500, 
                     n.ahead = 1, copula = "normal", alpha = 0.05, beta = 0.05)
t_VCoVaR = get_FC(data = cbind(rBTC, rLTC, rXMR, rXRP), type = "Vulnerability", window.length = 500, 
                  n.ahead = 1, copula = "t", alpha = 0.05, beta = 0.05)
Clayton_VCoVaR = get_FC(data = cbind(rBTC, rLTC, rXMR, rXRP), type = "Vulnerability", window.length = 500, 
                        n.ahead = 1, copula = "clayton", alpha = 0.05, beta = 0.05)
Gumbel_VCoVaR = get_FC(data = cbind(rBTC, rLTC, rXMR, rXRP), type = "Vulnerability", window.length = 500, 
                        n.ahead = 1, copula = "gumbel", alpha = 0.05, beta = 0.05)

###################################### EVALUATION #########################################

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

#VaR evaluation - Table 13
window.size = 500

length(which(rBTC[-c(1:window.size)] <= VaRBTC))/length(VaRBTC)
length(which(rLTC[-c(1:window.size)] <= VaRLTC))/length(VaRLTC)
length(which(rXMR[-c(1:window.size)] <= VaRXMR))/length(VaRXMR)
length(which(rXRP[-c(1:window.size)] <= VaRXRP))/length(VaRXRP)
length(which(rSys[-c(1:window.size)] <= VaRSys))/length(VaRSys)

#CoVaR evaluation - Table 14
norm_CoVaR_LTC_eval = evalBivariate(retX = rLTC[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                    VaRX = VaRLTC, CoVaR = norm_CoVaR_LTC)
norm_CoVaR_XMR_eval = evalBivariate(retX = rXMR[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                    VaRX = VaRXMR, CoVaR = norm_CoVaR_XMR)
norm_CoVaR_XRP_eval = evalBivariate(retX = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                    VaRX = VaRXRP, CoVaR = norm_CoVaR_XRP)
norm_SCoVaR_eval = evalBivariate(retX = rSys[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                 VaRX = VaRSys, CoVaR = norm_SCoVaR)

t_CoVaR_LTC_eval = evalBivariate(retX = rLTC[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                 VaRX = VaRLTC, CoVaR = t_CoVaR_LTC)
t_CoVaR_XMR_eval = evalBivariate(retX = rXMR[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                 VaRX = VaRXMR, CoVaR = t_CoVaR_XMR)
t_CoVaR_XRP_eval = evalBivariate(retX = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                 VaRX = VaRXRP, CoVaR = t_CoVaR_XRP)
t_SCoVaR_eval = evalBivariate(retX = rSys[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                              VaRX = VaRSys, CoVaR = t_SCoVaR)

Clayton_CoVaR_LTC_eval = evalBivariate(retX = rLTC[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                       VaRX = VaRLTC, CoVaR = Clayton_CoVaR_LTC)
Clayton_CoVaR_XMR_eval = evalBivariate(retX = rXMR[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                       VaRX = VaRXMR, CoVaR = Clayton_CoVaR_XMR)
Clayton_CoVaR_XRP_eval = evalBivariate(retX = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                       VaRX = VaRXRP, CoVaR = Clayton_CoVaR_XRP)
Clayton_SCoVaR_eval = evalBivariate(retX = rSys[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                    VaRX = VaRSys, CoVaR = Clayton_SCoVaR)

Gumbel_CoVaR_LTC_eval = evalBivariate(retX = rLTC[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                      VaRX = VaRLTC, CoVaR = Gumbel_CoVaR_LTC)
Gumbel_CoVaR_XMR_eval = evalBivariate(retX = rXMR[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                      VaRX = VaRXMR, CoVaR = Gumbel_CoVaR_XMR)
Gumbel_CoVaR_XRP_eval = evalBivariate(retX = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                      VaRX = VaRXRP, CoVaR = Gumbel_CoVaR_XRP)
Gumbel_SCoVaR_eval = evalBivariate(retX = rSys[-c(1:window.size)], retY = rBTC[-c(1:window.size)], 
                                   VaRX = VaRSys, CoVaR = Gumbel_SCoVaR)

norm_MCoVaR_eval = evalMulti3X(retX1 = rLTC[-c(1:window.size)], retX2 = rXMR[-c(1:window.size)],
                               retX3 = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)],
                               VaRX1 = VaRLTC, VaRX2 = VaRXMR, VaRX3 = VaRXRP, MCoVaR = norm_MCoVaR)
norm_VCoVaR_eval = evalVulnerability3X(retX1 = rLTC[-c(1:window.size)], retX2 = rXMR[-c(1:window.size)],
                                       retX3 = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)],
                                       VaRX1 = VaRLTC, VaRX2 = VaRXMR, VaRX3 = VaRXRP, VCoVaR = norm_VCoVaR)
 
t_MCoVaR_eval = evalMulti3X(retX1 = rLTC[-c(1:window.size)], retX2 = rXMR[-c(1:window.size)],
                            retX3 = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)],
                            VaRX1 = VaRLTC, VaRX2 = VaRXMR, VaRX3 = VaRXRP, MCoVaR = t_MCoVaR)
t_VCoVaR_eval = evalVulnerability3X(retX1 = rLTC[-c(1:window.size)], retX2 = rXMR[-c(1:window.size)],
                                    retX3 = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)],
                                    VaRX1 = VaRLTC, VaRX2 = VaRXMR, VaRX3 = VaRXRP, VCoVaR = t_VCoVaR)

Clayton_MCoVaR_eval = evalMulti3X(retX1 = rLTC[-c(1:window.size)], retX2 = rXMR[-c(1:window.size)],
                                  retX3 = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)],
                                  VaRX1 = VaRLTC, VaRX2 = VaRXMR, VaRX3 = VaRXRP, MCoVaR = Clayton_MCoVaR)
Clayton_VCoVaR_eval = evalVulnerability3X(retX1 = rLTC[-c(1:window.size)], retX2 = rXMR[-c(1:window.size)],
                                          retX3 = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)],
                                          VaRX1 = VaRLTC, VaRX2 = VaRXMR, VaRX3 = VaRXRP, VCoVaR = Clayton_VCoVaR)

Gumbel_MCoVaR_eval = evalMulti3X(retX1 = rLTC[-c(1:window.size)], retX2 = rXMR[-c(1:window.size)],
                                 retX3 = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)],
                                 VaRX1 = VaRLTC, VaRX2 = VaRXMR, VaRX3 = VaRXRP, MCoVaR = Gumbel_MCoVaR)
Gumbel_VCoVaR_eval = evalVulnerability3X(retX1 = rLTC[-c(1:window.size)], retX2 = rXMR[-c(1:window.size)],
                                         retX3 = rXRP[-c(1:window.size)], retY = rBTC[-c(1:window.size)],
                                         VaRX1 = VaRLTC, VaRX2 = VaRXMR, VaRX3 = VaRXRP, VCoVaR = Gumbel_VCoVaR)

res = rbind(c(norm_CoVaR_LTC_eval$rate, t_CoVaR_LTC_eval$rate, Clayton_CoVaR_LTC_eval$rate, Gumbel_CoVaR_LTC_eval$rate),
            c(norm_CoVaR_XMR_eval$rate, t_CoVaR_XMR_eval$rate, Clayton_CoVaR_XMR_eval$rate, Gumbel_CoVaR_XMR_eval$rate),
            c(norm_CoVaR_XRP_eval$rate, t_CoVaR_XRP_eval$rate, Clayton_CoVaR_XRP_eval$rate, Gumbel_CoVaR_XRP_eval$rate),
            c(norm_SCoVaR_eval$rate, t_SCoVaR_eval$rate, Clayton_SCoVaR_eval$rate, Gumbel_SCoVaR_eval$rate),
            c(norm_MCoVaR_eval$rate, t_MCoVaR_eval$rate, Clayton_MCoVaR_eval$rate, Gumbel_MCoVaR_eval$rate),
            c(norm_VCoVaR_eval$rate, t_VCoVaR_eval$rate, Clayton_VCoVaR_eval$rate, Gumbel_VCoVaR_eval$rate))
rownames(res) = c("CoVaR|LTC", "CoVaR|XMR", "CoVaR|XRP", "SCoVaR", "MCoVaR", "VCoVaR")
colnames(res) = c("normal", "t", "clayton", "gumbel")
