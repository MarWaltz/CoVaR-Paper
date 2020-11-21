#Note: The calculations in this script take approximately 15 minutes on an Intel(R) Xeon(R)
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

#Sys
specrSys = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,0)),
                      distribution.model = "sstd")
fitrSys = ugarchfit(spec = specrSys, data = rSys, solver = "hybrid")
u4 = as.vector(pit(fitrSys))


####################################### normal-COPULA ##############################################
#BTC-LTC
TIcop1 = fitCopula(copula = normalCopula(), data = cbind(v, u1), estimate.variance = F, method = "ml")
TIcop1@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrLTC)) + length(coef(TIcop1)))-2*(logLik(TIcop1) + fitrBTC@fit$LLH + fitrLTC@fit$LLH)) #AIC

#BTC-XMR
TIcop2 = fitCopula(copula = normalCopula(), data = cbind(v, u2), estimate.variance = F, method = "ml")
TIcop2@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrXMR)) + length(coef(TIcop2)))-2*(logLik(TIcop2) + fitrBTC@fit$LLH + fitrXMR@fit$LLH)) #AIC

#BTC-XRP
TIcop3 = fitCopula(copula = normalCopula(), data = cbind(v, u3), estimate.variance = F, method = "ml")
TIcop3@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrXRP)) + length(coef(TIcop3)))-2*(logLik(TIcop3) + fitrBTC@fit$LLH + fitrXRP@fit$LLH)) #AIC

#BTC-System
TIcop4 = fitCopula(copula = normalCopula(), data = cbind(v, u4), estimate.variance = F, method = "ml")
TIcop4@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrSys)) + length(coef(TIcop4)))-2*(logLik(TIcop4) + fitrBTC@fit$LLH + fitrSys@fit$LLH)) #AIC

#BTC-LTC-XMR-XRP
TIfitMulti = fitCopula(copula = normalCopula(dim = 4, dispstr = "un"), estimate.variance = F,
                        data = cbind(v, u1, u2, u3), method = "ml")
p2P(TIfitMulti@estimate) #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrLTC)) + length(coef(fitrXMR)) + length(coef(fitrXRP)) + length(coef(TIfitMulti)))
           -2*(logLik(TIfitMulti) + fitrBTC@fit$LLH + fitrLTC@fit$LLH + fitrXMR@fit$LLH + fitrXRP@fit$LLH)) #AIC


########################################### t-COPULA ###############################################
#BTC-LTC
TIcop1 = fitCopula(copula = tCopula(), data = cbind(v, u1), estimate.variance = F, method = "ml")
TIcop1@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrLTC)) + length(coef(TIcop1)))-2*(logLik(TIcop1) + fitrBTC@fit$LLH + fitrLTC@fit$LLH)) #AIC

#BTC-XMR
TIcop2 = fitCopula(copula = tCopula(), data = cbind(v, u2), estimate.variance = F, method = "ml")
TIcop2@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrXMR)) + length(coef(TIcop2)))-2*(logLik(TIcop2) + fitrBTC@fit$LLH + fitrXMR@fit$LLH)) #AIC

#BTC-XRP
TIcop3 = fitCopula(copula = tCopula(), data = cbind(v, u3), estimate.variance = F, method = "ml")
TIcop3@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrXRP)) + length(coef(TIcop3)))-2*(logLik(TIcop3) + fitrBTC@fit$LLH + fitrXRP@fit$LLH)) #AIC

#BTC-System
TIcop4 = fitCopula(copula = tCopula(), data = cbind(v, u4), estimate.variance = F, method = "ml")
TIcop4@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrSys)) + length(coef(TIcop4)))-2*(logLik(TIcop4) + fitrBTC@fit$LLH + fitrSys@fit$LLH)) #AIC

#BTC-LTC-XMR-XRP
TIfitMulti = fitCopula(copula = tCopula(dim = 4, dispstr = "un"), estimate.variance = F,
                        data = cbind(v, u1, u2, u3), method = "ml")
p2P(TIfitMulti@estimate[-7]) #copula estimate
TIfitMulti@estimate[7] #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrLTC)) + length(coef(fitrXMR)) + length(coef(fitrXRP)) + length(coef(TIfitMulti)))
           -2*(logLik(TIfitMulti) + fitrBTC@fit$LLH + fitrLTC@fit$LLH + fitrXMR@fit$LLH + fitrXRP@fit$LLH)) #AIC

######################################## clayton-COPULA #############################################
#BTC-LTC
TIcop1 = fitCopula(copula = claytonCopula(), data = cbind(v, u1), estimate.variance = F, method = "ml")
TIcop1@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrLTC)) + length(coef(TIcop1)))-2*(logLik(TIcop1) + fitrBTC@fit$LLH + fitrLTC@fit$LLH)) #AIC

#BTC-XMR
TIcop2 = fitCopula(copula = claytonCopula(), data = cbind(v, u2), estimate.variance = F, method = "ml")
TIcop2@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrXMR)) + length(coef(TIcop2)))-2*(logLik(TIcop2) + fitrBTC@fit$LLH + fitrXMR@fit$LLH)) #AIC

#BTC-XRP
TIcop3 = fitCopula(copula = claytonCopula(), data = cbind(v, u3), estimate.variance = F, method = "itau") #ML estimation failes
TIcop3@estimate #copula estimate
#AIC is not computed as Maximum Likelihood fails

#BTC-System
TIcop4 = fitCopula(copula = claytonCopula(), data = cbind(v, u4), estimate.variance = F, method = "ml")
TIcop4@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrSys)) + length(coef(TIcop4)))-2*(logLik(TIcop4) + fitrBTC@fit$LLH + fitrSys@fit$LLH)) #AIC

#BTC-LTC-XMR-XRP
TIfitMulti = fitCopula(copula = claytonCopula(dim = 4), estimate.variance = F,
                        data = cbind(v, u1, u2, u3), method = "ml")
TIfitMulti@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrLTC)) + length(coef(fitrXMR)) + length(coef(fitrXRP)) + length(coef(TIfitMulti)))
           -2*(logLik(TIfitMulti) + fitrBTC@fit$LLH + fitrLTC@fit$LLH + fitrXMR@fit$LLH + fitrXRP@fit$LLH)) #AIC

######################################## gumbel-COPULA #############################################
#BTC-LTC
TIcop1 = fitCopula(copula = gumbelCopula(), data = cbind(v, u1), estimate.variance = F, method = "ml")
TIcop1@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrLTC)) + length(coef(TIcop1)))-2*(logLik(TIcop1) + fitrBTC@fit$LLH + fitrLTC@fit$LLH)) #AIC

#BTC-XMR
TIcop2 = fitCopula(copula = gumbelCopula(), data = cbind(v, u2), estimate.variance = F, method = "ml")
TIcop2@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrXMR)) + length(coef(TIcop2)))-2*(logLik(TIcop2) + fitrBTC@fit$LLH + fitrXMR@fit$LLH)) #AIC

#BTC-XRP
TIcop3 = fitCopula(copula = gumbelCopula(), data = cbind(v, u3), estimate.variance = F, method = "ml")
TIcop3@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrXRP)) + length(coef(TIcop3)))-2*(logLik(TIcop3) + fitrBTC@fit$LLH + fitrXRP@fit$LLH)) #AIC

#BTC-System
TIcop4 = fitCopula(copula = gumbelCopula(), data = cbind(v, u4), estimate.variance = F, method = "ml")
TIcop4@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrSys)) + length(coef(TIcop4)))-2*(logLik(TIcop4) + fitrBTC@fit$LLH + fitrSys@fit$LLH)) #AIC

#BTC-LTC-XMR-XRP
TIfitMulti = fitCopula(copula = gumbelCopula(dim = 4), estimate.variance = F,
                        data = cbind(v, u1, u2, u3), method = "ml")
TIfitMulti@estimate #copula estimate
as.numeric(2*(length(coef(fitrBTC)) + length(coef(fitrLTC)) + length(coef(fitrXMR)) + length(coef(fitrXRP)) + length(coef(TIfitMulti)))
           -2*(logLik(TIfitMulti) + fitrBTC@fit$LLH + fitrLTC@fit$LLH + fitrXMR@fit$LLH + fitrXRP@fit$LLH)) #AIC

####################################### Patton-t-COPULA ############################################
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

### Optimization
#BTC-LTC
out_BTC_LTC = nloptr(x0 = c(0.5, 0.5, 0.5, 4), eval_f = ObjectiveF, lb = c(-Inf, -Inf, -Inf, 2),
                     ub = c(Inf, Inf, Inf, 200), opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-04),
                     data = cbind(v, u1))
round(out_BTC_LTC$solution,4) #parameter estimates in order: omega_theta, beta_theta, c_theta, v
2*(length(coef(fitrBTC)) + length(coef(fitrLTC)) + length(out_BTC_LTC$solution))-2*(-out_BTC_LTC$objective + fitrBTC@fit$LLH + fitrLTC@fit$LLH) #AIC

#BTC-XMR
out_BTC_XMR = nloptr(x0 = c(0.5, 0.5, 0.5, 4), eval_f = ObjectiveF, lb = c(-Inf, -Inf, -Inf, 2),
                     ub = c(Inf, Inf, Inf, 200), opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-04),
                     data = cbind(v, u2))
round(out_BTC_XMR$solution,4) #parameter estimates in order: omega_theta, beta_theta, c_theta, v
2*(length(coef(fitrBTC)) + length(coef(fitrXMR)) + length(out_BTC_XMR$solution))-2*(-out_BTC_XMR$objective + fitrBTC@fit$LLH + fitrXMR@fit$LLH) #AIC

#BTC-XRP
out_BTC_XRP = nloptr(x0 = c(0.5, 0.5, 0.5, 4), eval_f = ObjectiveF, lb = c(-Inf, -Inf, -Inf, 2),
                     ub = c(Inf, Inf, Inf, 200), opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-04),
                     data = cbind(v, u3))
round(out_BTC_XRP$solution,4) #parameter estimates in order: omega_theta, beta_theta, c_theta, v
2*(length(coef(fitrBTC)) + length(coef(fitrXRP)) + length(out_BTC_XRP$solution))-2*(-out_BTC_XRP$objective + fitrBTC@fit$LLH + fitrXRP@fit$LLH) #AIC

#BTC-System
out_BTC_Sys = nloptr(x0 = c(0.5, 0.5, 0.5, 4), eval_f = ObjectiveF, lb = c(-Inf, -Inf, -Inf, 2),
                     ub = c(Inf, Inf, Inf, 200), opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-04),
                     data = cbind(v, u4))
round(out_BTC_Sys$solution,4) #parameter estimates in order: omega_theta, beta_theta, c_theta, v
2*(length(coef(fitrBTC)) + length(coef(fitrSys)) + length(out_BTC_Sys$solution))-2*(-out_BTC_Sys$objective + fitrBTC@fit$LLH + fitrSys@fit$LLH) #AIC


######################################### DCC-t-COPULA ############################################
#BTC-LTC
BivDCCspec1 = cgarchspec(uspec = multispec(list(specrBTC, specrLTC)), dccOrder = c(1,1),
                         distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                   transformation = "parametric"))
BivDCCfit1 = cgarchfit(spec = BivDCCspec1, data = cbind(rBTC, rLTC), solver = c("hybrid", "solnp"))
coef(BivDCCfit1, type = "dcc") #copula estimates
2*length(coef(BivDCCfit1)) -2*(likelihood(BivDCCfit1)) #AIC

#BTC-XMR
BivDCCspec2 = cgarchspec(uspec = multispec(list(specrBTC, specrXMR)), dccOrder = c(1,1),
                         distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                   transformation = "parametric"))
BivDCCfit2 = cgarchfit(spec = BivDCCspec2, data = cbind(rBTC, rXMR), solver = c("hybrid", "solnp"))
coef(BivDCCfit2, type = "dcc") #copula estimates
2*length(coef(BivDCCfit2)) -2*(likelihood(BivDCCfit2)) #AIC

#BTC-XRP
BivDCCspec3 = cgarchspec(uspec = multispec(list(specrBTC, specrXRP)), dccOrder = c(1,1),
                         distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                   transformation = "parametric"))
BivDCCfit3 = cgarchfit(spec = BivDCCspec3, data = cbind(rBTC, rXRP), solver = c("hybrid", "solnp"))
coef(BivDCCfit3, type = "dcc") #copula estimates
2*length(coef(BivDCCfit3)) -2*(likelihood(BivDCCfit3)) #AIC

#BTC-System
BivDCCspec4 = cgarchspec(uspec = multispec(list(specrBTC, specrSys)), dccOrder = c(1,1),
                         distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                   transformation = "parametric"))
BivDCCfit4 = cgarchfit(spec = BivDCCspec4, data = cbind(rBTC, rSys), solver = c("hybrid", "solnp"))
coef(BivDCCfit4, type = "dcc") #copula estimates
2*length(coef(BivDCCfit4)) -2*(likelihood(BivDCCfit4)) #AIC

#BTC-LTC-XMR-XRP
DCCspec1 = cgarchspec(uspec = multispec(list(specrBTC, specrLTC, specrXMR, specrXRP)), dccOrder = c(1,1),
                      distribution.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                transformation = "parametric"))
DCCfit1 = cgarchfit(spec = DCCspec1, data = cbind(rBTC, rLTC, rXMR, rXRP), solver = c("hybrid", "solnp"))
coef(DCCfit1, type = "dcc") #copula estimates
2*length(coef(DCCfit1)) -2*(likelihood(DCCfit1)) #AIC
