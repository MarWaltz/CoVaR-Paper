#Note: The calculations in this script takes only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")

#load BTC, LTC, XMR, XRP from 'https://github.com/MarWaltz/CoVaR-Paper/tree/master/Section 5/Data'
rBTC = diff(log(BTC))
rLTC = diff(log(LTC))
rXMR = diff(log(XMR))
rXRP = diff(log(XRP))
rSys = rLTC + rXMR + rXRP

set.seed(312)

#################################### MODEL SPECIFICATION ##########################################
#BTC
specrBTC = ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean = 1),
                      variance.model = list(model = "gjrGARCH", garchOrder = c(3,0)),
                      distribution.model = "sstd")
fitrBTC = ugarchfit(spec = specrBTC, data = rBTC, solver = "hybrid")
resBTC = as.vector(residuals(fitrBTC, standardize = T))

#LTC
specrLTC = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrLTC = ugarchfit(spec = specrLTC, data = rLTC, solver = "hybrid")
resLTC = as.vector(residuals(fitrLTC, standardize = T))

#XMR
specrXMR = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,1)),
                      distribution.model = "sstd")
fitrXMR = ugarchfit(spec = specrXMR, data = rXMR, solver = "hybrid")
resXMR = as.vector(residuals(fitrXMR, standardize = T))

#XRP
specrXRP = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 1),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrXRP = ugarchfit(spec = specrXRP, data = rXRP, solver = "hybrid")
resXRP = as.vector(residuals(fitrXRP, standardize = T))

#System
specrSys = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,0)),
                      distribution.model = "sstd")
fitrSys = ugarchfit(spec = specrSys, data = rSys, solver = "hybrid")
resSys = as.vector(residuals(fitrSys, standardize = T))


####################################### ACF/PACF PLOTS ############################################
acf(resBTC, main = "ACF - BTC")
pacf(resBTC, main = "PACF - BTC")

acf(resLTC, main = "ACF - LTC")
pacf(resLTC, main = "PACF - LTC")

acf(resXMR, main = "ACF - XMR")
pacf(resXMR, main = "PACF - XMR")

acf(resXRP, main = "ACF - XRP")
pacf(resXRP, main = "PACF - XRP")

acf(resSys, main = "ACF - System")
pacf(resSys, main = "PACF - System")
