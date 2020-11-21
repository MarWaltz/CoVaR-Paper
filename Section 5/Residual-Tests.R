#Note: The calculations in this script take only a couple of seconds/minutes on an Intel(R) 
#      Xeon(R) Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")
library("WeightedPortTest")

#load BTC, LTC, XMR, XRP from 'https://github.com/MarWaltz/Thesis/tree/master/Section 5/Data'
rBTC = diff(log(BTC))
rLTC = diff(log(LTC))
rXMR = diff(log(XMR))
rXRP = diff(log(XRP))
rSys = rLTC + rXMR + rXRP

set.seed(312)

#BTC
specrBTC = ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean = 1),
                      variance.model = list(model = "gjrGARCH", garchOrder = c(3,0)),
                      distribution.model = "sstd")
fitrBTC = ugarchfit(spec = specrBTC, data = rBTC, solver = "hybrid")
colBTC = round(c("AIC" = infocriteria(fitrBTC)[1]*length(rBTC), #'rugarch' computes AIC per observation
                 "LB" = Box.test(as.numeric(residuals(fitrBTC, standardize = T)), lag = 8, 
                                 fitdf = 2, type = c("Ljung-Box"))$p.value,
                 "LBsq" = Box.test(as.numeric(residuals(fitrBTC, standardize = T))^2, lag = 8, 
                                   fitdf = 3, type = c("Ljung-Box"))$p.value,
                 "WLM" = Weighted.LM.test(as.numeric(residuals(fitrBTC)), h.t = as.numeric(sigma(fitrBTC))^2, 
                                          weighted = T, lag = 8, type = "correlation", fitdf = 3)$p.value,
                 "Sign Bias" = signbias(fitrBTC)[1,2],
                 "Neg. Sign Bias" = signbias(fitrBTC)[2,2],
                 "Pos. Sign Bias" = signbias(fitrBTC)[3,2],
                 "Joint effect"= signbias(fitrBTC)[4,2]),4)

#LTC
specrLTC = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrLTC = ugarchfit(spec = specrLTC, data = rLTC, solver = "hybrid")
colLTC = round(c("AIC" = infocriteria(fitrLTC)[1]*length(rLTC),
                 "LB" = Box.test(as.numeric(residuals(fitrLTC, standardize = T)), lag = 8, 
                                 fitdf = 0, type = c("Ljung-Box"))$p.value,
                 "LBsq" = Box.test(as.numeric(residuals(fitrLTC, standardize = T))^2, lag = 8, 
                                   fitdf = 2, type = c("Ljung-Box"))$p.value,
                 "WLM" = Weighted.LM.test(as.numeric(residuals(fitrLTC)), h.t = as.numeric(sigma(fitrLTC))^2, 
                                          weighted = T, lag = 8, type = "correlation", fitdf = 2)$p.value,
                 "Sign Bias" = signbias(fitrLTC)[1,2],
                 "Neg. Sign Bias" = signbias(fitrLTC)[2,2],
                 "Pos. Sign Bias" = signbias(fitrLTC)[3,2],
                 "Joint effect"= signbias(fitrLTC)[4,2]),4)

#XMR
specrXMR = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,1)),
                      distribution.model = "sstd")
fitrXMR = ugarchfit(spec = specrXMR, data = rXMR, solver = "hybrid")
colXMR = round(c("AIC" = infocriteria(fitrXMR)[1]*length(rXMR),
                 "LB" = Box.test(as.numeric(residuals(fitrXMR, standardize = T)), lag = 8, 
                                 fitdf = 0, type = c("Ljung-Box"))$p.value,
                 "LBsq" = Box.test(as.numeric(residuals(fitrXMR, standardize = T))^2, lag = 8, 
                                   fitdf = 7, type = c("Ljung-Box"))$p.value,
                 "WLM" = Weighted.LM.test(as.numeric(residuals(fitrXMR)), h.t = as.numeric(sigma(fitrXMR))^2, 
                                          weighted = T, lag = 8, type = "correlation", fitdf = 7)$p.value,
                 "Sign Bias" = signbias(fitrXMR)[1,2],
                 "Neg. Sign Bias" = signbias(fitrXMR)[2,2],
                 "Pos. Sign Bias" = signbias(fitrXMR)[3,2],
                 "Joint effect"= signbias(fitrXMR)[4,2]),4)

#XRP
specrXRP = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 1),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrXRP = ugarchfit(spec = specrXRP, data = rXRP, solver = "hybrid")
colXRP = round(c("AIC" = infocriteria(fitrXRP)[1]*length(rXRP),
                 "LB" = Box.test(as.numeric(residuals(fitrXRP, standardize = T)), lag = 8, 
                                 fitdf = 0, type = c("Ljung-Box"))$p.value,
                 "LBsq" = Box.test(as.numeric(residuals(fitrXRP, standardize = T))^2, lag = 8, 
                                   fitdf = 2, type = c("Ljung-Box"))$p.value,
                 "WLM" = Weighted.LM.test(as.numeric(residuals(fitrXRP)), h.t = as.numeric(sigma(fitrXRP))^2, 
                                          weighted = T, lag = 8, type = "correlation", fitdf = 2)$p.value,
                 "Sign Bias" = signbias(fitrXRP)[1,2],
                 "Neg. Sign Bias" = signbias(fitrXRP)[2,2],
                 "Pos. Sign Bias" = signbias(fitrXRP)[3,2],
                 "Joint effect"= signbias(fitrXRP)[4,2]),4)

#System
specrSys = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,0)),
                      distribution.model = "sstd")
fitrSys = ugarchfit(spec = specrSys, data = rSys, solver = "hybrid")
colSys = round(c("AIC" = infocriteria(fitrSys)[1]*length(rSys),
                 "LB" = Box.test(as.numeric(residuals(fitrSys, standardize = T)), lag = 8, 
                                 fitdf = 0, type = c("Ljung-Box"))$p.value,
                 "LBsq" = Box.test(as.numeric(residuals(fitrSys, standardize = T))^2, lag = 8, 
                                   fitdf = 6, type = c("Ljung-Box"))$p.value,
                 "WLM" = Weighted.LM.test(as.numeric(residuals(fitrSys)), h.t = as.numeric(sigma(fitrSys))^2, 
                                          weighted = T, lag = 8, type = "correlation", fitdf = 6)$p.value,
                 "Sign Bias" = signbias(fitrSys)[1,2],
                 "Neg. Sign Bias" = signbias(fitrSys)[2,2],
                 "Pos. Sign Bias" = signbias(fitrSys)[3,2],
                 "Joint effect"= signbias(fitrSys)[4,2]),4)

res = cbind("BTC" = colBTC, "LTC" = colLTC, "XMR" = colXMR, "XRP" = colXRP, "System" = colSys)
