#Note: The calculations in this script take approximately 7h on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")
library("WeightedPortTest")

#load BTC, LTC, XMR, XRP from 'https://github.com/MarWaltz/Thesis/tree/master/Section 5/Data'
rBTC = diff(log(BTC))
rLTC = diff(log(LTC))
rXMR = diff(log(XMR))
rXRP = diff(log(XRP))
rSys = rLTC + rXMR + rXRP

getAllARMA = function(ret, ar.max = 6, ma.max = 6){
  table = expand.grid(0:ar.max, 0:ma.max, 0:1)
  extraCols = 3
  for(i in seq_len(extraCols)){
    table = cbind(table, rep(NA, nrow(table)))  
  }
  colnames(table) = c("AR", "MA", "Mean", "AIC", "LB", "LBsq")
  
  for(i in seq_len(nrow(table))){
    #spec and fit
    spec = arfimaspec(mean.model = list(armaOrder = c(table$AR[i], table$MA[i]), include.mean = table$Mean[i]),
                      distribution.model = "sstd")
    
    fit = try(arfimafit(spec, data = as.vector(ret), solver = "hybrid"), silent = T)
    if(class(fit) == "try-error"){next}
    
    #AIC
    table$AIC[i] = infocriteria(fit)[1]
    
    #tests
    res = as.vector(residuals(fit))
    df = sum(c(table$AR[i], table$MA[i]))
    
    tmp = try(Box.test(res, lag = max(df+1, 8), fitdf = df, type = c("Ljung-Box"))$p.value, silent = T)
    if(class(tmp)!= "try-error"){
      table$LB[i] = tmp  
    }
    
    tmp = try(Box.test(res^2, lag = max(df+1, 8), fitdf = df, type = c("Ljung-Box"))$p.value, silent = T)
    if(class(tmp)!= "try-error"){
      table$LBsq[i] = tmp  
    }
    
    res = df = NULL
    cat(paste("Iteration done: ", i, "/", nrow(table), "\n"))
  }
  return(table)
}

getAllARMAGARCH = function(ret, ar.max = 6, ma.max = 6, arch.max = 6, garch.max = 1, model = "sGARCH"){
  table = expand.grid(0:ar.max, 0:ma.max, 0:1, 0:arch.max, 0:garch.max)
  extraCols = 3
  for(i in seq_len(extraCols)){
    table = cbind(table, rep(NA, nrow(table)))  
  }
  colnames(table) = c("AR", "MA", "Mean", "ARCH", "GARCH", "AIC", "LB", "WLM")
  
  for(i in seq_len(nrow(table))){
    #spec and fit
    if(model == "sGARCH"){
      spec = ugarchspec(mean.model = list(armaOrder = c(table$AR[i], table$MA[i]),include.mean = table$Mean[i]),
                        variance.model = list(model = "sGARCH", garchOrder = c(table$ARCH[i], table$GARCH[i])),
                        distribution.model = "sstd")
    }
    if(model == "gjrGARCH"){
      spec = ugarchspec(mean.model = list(armaOrder = c(table$AR[i], table$MA[i]),include.mean = table$Mean[i]),
                        variance.model = list(model = "gjrGARCH", garchOrder = c(table$ARCH[i], table$GARCH[i])),
                        distribution.model = "sstd")     
    }
    
    fit = try(ugarchfit(spec, data = as.vector(ret), solver = "hybrid"), silent = T)
    if(class(fit) == "try-error"){next}
    
    #AIC
    table$AIC[i] = infocriteria(fit)[1]
    
    #tests
    res = as.numeric(residuals(fit))
    stdres = as.numeric(residuals(fit, standardize = TRUE))
    h.t = as.numeric(sigma(fit))^2
    df = sum(c(table$AR[i], table$MA[i]))
    gdf = sum(c(table$ARCH[i], table$GARCH[i]))
    
    tmp = try(Box.test(stdres, lag = max(df+1, 8), fitdf = df, type = c("Ljung-Box"))$p.value, silent = T)
    if(class(tmp)!= "try-error"){
      table$LB[i] = tmp  
    }
    
    tmp = try(Weighted.LM.test(res, h.t = h.t, weighted = T, lag = max(gdf+1, 8), type = "correlation", fitdf = gdf)$p.value, silent = T)
    if(class(tmp)!= "try-error"){
      table$WLM[i] = tmp  
    }
    
    res = stdres = h.t = df = gdf = NULL
    cat(paste("Iteration done: ", i, "/", nrow(table), "\n"))
  }
  return(table)
}

set.seed(127)

############# Step 1: Calculate all ARMA(p,q) models and consider those fulfilling the LB #########
BTC_ARMA = getAllARMA(rBTC, ar.max = 6, ma.max = 6)
BTC_ARMA_red = BTC_ARMA[which(BTC_ARMA$LB > 0.05),]

LTC_ARMA = getAllARMA(rLTC, ar.max = 6, ma.max = 6)
LTC_ARMA_red = LTC_ARMA[which(LTC_ARMA$LB > 0.05),]

XMR_ARMA = getAllARMA(rXMR, ar.max = 6, ma.max = 6)
XMR_ARMA_red = XMR_ARMA[which(XMR_ARMA$LB > 0.05),]

XRP_ARMA = getAllARMA(rXRP, ar.max = 6, ma.max = 6)
XRP_ARMA_red = XRP_ARMA[which(XRP_ARMA$LB > 0.05),]

SYS_ARMA = getAllARMA(rSys, ar.max = 6, ma.max = 6)
SYS_ARMA_red = SYS_ARMA[which(SYS_ARMA$LB > 0.05),]

################# Step 2: Select best model according to AIC and test for ARCH effects ############
BTC_OPT = BTC_ARMA_red[order(BTC_ARMA_red$AIC),][1,] #0-1-1
#ARCH effects existent: Go to Step 3

LTC_OPT = LTC_ARMA_red[order(LTC_ARMA_red$AIC),][1,] 
#No LB test successful: Go to Step 3

XMR_OPT = XMR_ARMA_red[order(XMR_ARMA_red$AIC),][1,]
#No LB test successful: Go to Step 3

XRP_OPT = XRP_ARMA_red[order(XRP_ARMA_red$AIC),][1,]
#No LB test successful: Go to Step 3

SYS_OPT = SYS_ARMA_red[order(SYS_ARMA_red$AIC),][1,]
#No LB test successful: Go to Step 3


### Step 3: Calculate all ARMA(p,q)-GARCH(P,Q) models and consider those fulfilling LB, WLM ####
BTC_ARMA_GARCH = getAllARMAGARCH(rBTC, ar.max = 6, ma.max = 6, arch.max = 6, garch.max = 1, model = "sGARCH")
BTC_ARMA_GARCH_red = BTC_ARMA_GARCH[which((BTC_ARMA_GARCH$LB > 0.05) & (BTC_ARMA_GARCH$WLM > 0.05)),]

LTC_ARMA_GARCH = getAllARMAGARCH(rLTC, ar.max = 6, ma.max = 6, arch.max = 6, garch.max = 1, model = "sGARCH")
LTC_ARMA_GARCH_red = LTC_ARMA_GARCH[which((LTC_ARMA_GARCH$LB > 0.05) & (LTC_ARMA_GARCH$WLM > 0.05)),]

XMR_ARMA_GARCH = getAllARMAGARCH(rXMR, ar.max = 6, ma.max = 6, arch.max = 6, garch.max = 1, model = "sGARCH")
XMR_ARMA_GARCH_red = XMR_ARMA_GARCH[which((XMR_ARMA_GARCH$LB > 0.05) & (XMR_ARMA_GARCH$WLM > 0.05)),]

XRP_ARMA_GARCH = getAllARMAGARCH(rXRP, ar.max = 6, ma.max = 6, arch.max = 6, garch.max = 1, model = "sGARCH")
XRP_ARMA_GARCH_red = XRP_ARMA_GARCH[which((XRP_ARMA_GARCH$LB > 0.05) & (XRP_ARMA_GARCH$WLM > 0.05)),]

SYS_ARMA_GARCH = getAllARMAGARCH(rSys, ar.max = 6, ma.max = 6, arch.max = 6, garch.max = 1, model = "sGARCH")
SYS_ARMA_GARCH_red = SYS_ARMA_GARCH[which((SYS_ARMA_GARCH$LB > 0.05) & (SYS_ARMA_GARCH$WLM > 0.05)),]

############# Step 4: Select best model according to AIC and test for leverage effects #############
BTC_OPT = BTC_ARMA_GARCH_red[order(BTC_ARMA_GARCH_red$AIC),][1,]

LTC_OPT = LTC_ARMA_GARCH_red[order(LTC_ARMA_GARCH_red$AIC),][1,]

XMR_OPT = XMR_ARMA_GARCH_red[order(XMR_ARMA_GARCH_red$AIC),][1,]

XRP_OPT = XRP_ARMA_GARCH_red[order(XRP_ARMA_GARCH_red$AIC),][1,]

SYS_OPT = SYS_ARMA_GARCH_red[order(SYS_ARMA_GARCH_red$AIC),][1,]

#BTC
specY = ugarchspec(mean.model = list(armaOrder = as.numeric(BTC_OPT[1:2]),include.mean = as.numeric(BTC_OPT[3])),
                   variance.model = list(model = "sGARCH", garchOrder = as.numeric(BTC_OPT[4:5])),
                   distribution.model = "sstd")
fitBTC = ugarchfit(spec = specY, data = rBTC, solver = "hybrid")
signbias(fitBTC) 
#fails at 5% level: Go to Step 5

#LTC
specY = ugarchspec(mean.model = list(armaOrder = as.numeric(LTC_OPT[1:2]),include.mean = as.numeric(LTC_OPT[3])),
                   variance.model = list(model = "sGARCH", garchOrder = as.numeric(LTC_OPT[4:5])),
                   distribution.model = "sstd")
fitLTC = ugarchfit(spec = specY, data = rLTC, solver = "hybrid")
signbias(fitLTC)
#fine: Final model is found

#XMR
specY = ugarchspec(mean.model = list(armaOrder = as.numeric(XMR_OPT[1:2]),include.mean = as.numeric(XMR_OPT[3])),
                   variance.model = list(model = "sGARCH", garchOrder = as.numeric(XMR_OPT[4:5])),
                   distribution.model = "sstd")
fitXMR = ugarchfit(spec = specY, data = rXMR, solver = "hybrid")
signbias(fitXMR)
#fine: Final model is found

#XRP
specY = ugarchspec(mean.model = list(armaOrder = as.numeric(XRP_OPT[1:2]),include.mean = as.numeric(XRP_OPT[3])),
                   variance.model = list(model = "sGARCH", garchOrder = as.numeric(XRP_OPT[4:5])),
                   distribution.model = "sstd")
fitXRP = ugarchfit(spec = specY, data = rXRP, solver = "hybrid")
signbias(fitXRP)
#fine: Final model is found

#Sys
specY = ugarchspec(mean.model = list(armaOrder = as.numeric(SYS_OPT[1:2]),include.mean = as.numeric(SYS_OPT[3])),
                   variance.model = list(model = "sGARCH", garchOrder = as.numeric(SYS_OPT[4:5])),
                   distribution.model = "sstd")
fitSys = ugarchfit(spec = specY, data = rSys, solver = "hybrid")
signbias(fitSys)
#fine: Final model is found

## Steps 5 & 6: Calculate all possible ARMA(p,q)-GJR-GARCH(P,Q) models, select the one with best AIC fulfilling LB, WLM, SB ##
BTC_ARMA_GJR_GARCH = getAllARMAGARCH(rBTC, ar.max = 6, ma.max = 6, arch.max = 6, garch.max = 1, model = "gjrGARCH")
BTC_ARMA_GJR_GARCH_red = BTC_ARMA_GJR_GARCH[which((BTC_ARMA_GJR_GARCH$LB > 0.05) & (BTC_ARMA_GJR_GARCH$WLM > 0.05)),]
BTC_OPT = BTC_ARMA_GJR_GARCH_red[order(BTC_ARMA_GJR_GARCH_red$AIC),][1,]

#BTC
specY = ugarchspec(mean.model = list(armaOrder = as.numeric(BTC_OPT[1:2]),include.mean = as.numeric(BTC_OPT[3])),
                   variance.model = list(model = "gjrGARCH", garchOrder = as.numeric(BTC_OPT[4:5])),
                   distribution.model = "sstd")
fitBTC = ugarchfit(spec = specY, data = rBTC, solver = "hybrid")
signbias(fitBTC) 
#fine: Final model is found


#Optimal Models:
#BTC: ARMA(1,1)-GJR-GARCH(3,0) with mu
#LTC: GARCH(1,1) without mu
#XMR: GARCH(6,1) without mu
#XRP: GARCH(1,1) with mu
#SYS: GARCH(6,0) without mu

#Parameters
coef(fitBTC)
coef(fitLTC)
coef(fitXMR)
coef(fitXRP)
coef(fitSys)

#Standard errors
round(fitBTC@fit$se.coef,4)
round(fitLTC@fit$se.coef,4)
round(fitXMR@fit$se.coef,4)
round(fitXRP@fit$se.coef,4)
round(fitSys@fit$se.coef,4)
