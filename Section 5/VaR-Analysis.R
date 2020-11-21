#Note: The calculations in this script take only a couple of seconds/minutes on an Intel(R) 
#      Xeon(R) Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")

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

#LTC
specrLTC = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrLTC = ugarchfit(spec = specrLTC, data = rLTC, solver = "hybrid")

#XMR
specrXMR = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,1)),
                      distribution.model = "sstd")
fitrXMR = ugarchfit(spec = specrXMR, data = rXMR, solver = "hybrid")

#XRP
specrXRP = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 1),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrXRP = ugarchfit(spec = specrXRP, data = rXRP, solver = "hybrid")

#System
specrSys = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,0)),
                      distribution.model = "sstd")
fitrSys = ugarchfit(spec = specrSys, data = rSys, solver = "hybrid")

####################################### Value at Risk ############################################
VaR_rBTC = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", 0.05, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])
VaR_rLTC = as.vector(fitted(fitrLTC)) + as.vector(sigma(fitrLTC))*qdist("sstd", 0.05, shape = coef(fitrLTC)["shape"], skew = coef(fitrLTC)["skew"])
VaR_rXMR = as.vector(fitted(fitrXMR)) + as.vector(sigma(fitrXMR))*qdist("sstd", 0.05, shape = coef(fitrXMR)["shape"], skew = coef(fitrXMR)["skew"])
VaR_rXRP = as.vector(fitted(fitrXRP)) + as.vector(sigma(fitrXRP))*qdist("sstd", 0.05, shape = coef(fitrXRP)["shape"], skew = coef(fitrXRP)["skew"])
VaR_rSys = as.vector(fitted(fitrSys)) + as.vector(sigma(fitrSys))*qdist("sstd", 0.05, shape = coef(fitrSys)["shape"], skew = coef(fitrSys)["skew"])

#test via internal VaR calculation of 'rugarch': 
#all(VaR_rBTC == quantile(fitrBTC, 0.05)) #TRUE

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

#aggregate results
res = cbind(rep(0.05,5),
            round(c(VaR_BTC_exc_rate, VaR_LTC_exc_rate, VaR_XMR_exc_rate,
                    VaR_XRP_exc_rate, VaR_Sys_exc_rate),4),
            rep(length(rBTC),5)*0.05,
            c(length(VaR_BTC_exc), length(VaR_LTC_exc), length(VaR_XMR_exc), 
              length(VaR_XRP_exc), length(VaR_Sys_exc)))
rownames(res) = c("BTC", "LTC", "XMR", "XRP", "System")
colnames(res) = c("alpha-level", "Measured Rate", "Theoretical Violations", "Measured Violations")
