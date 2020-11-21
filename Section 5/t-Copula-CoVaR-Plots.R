#Note: The calculations in this script take only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")
library("copula")
library("zoo")
library("lubridate")

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

#Univariate VaRs
VaR_rBTC = as.vector(fitted(fitrBTC)) + as.vector(sigma(fitrBTC))*qdist("sstd", 0.05, shape = coef(fitrBTC)["shape"], skew = coef(fitrBTC)["skew"])
VaR_rLTC = as.vector(fitted(fitrLTC)) + as.vector(sigma(fitrLTC))*qdist("sstd", 0.05, shape = coef(fitrLTC)["shape"], skew = coef(fitrLTC)["skew"])
VaR_rXMR = as.vector(fitted(fitrXMR)) + as.vector(sigma(fitrXMR))*qdist("sstd", 0.05, shape = coef(fitrXMR)["shape"], skew = coef(fitrXMR)["skew"])
VaR_rXRP = as.vector(fitted(fitrXRP)) + as.vector(sigma(fitrXRP))*qdist("sstd", 0.05, shape = coef(fitrXRP)["shape"], skew = coef(fitrXRP)["skew"])
VaR_rSys = as.vector(fitted(fitrSys)) + as.vector(sigma(fitrSys))*qdist("sstd", 0.05, shape = coef(fitrSys)["shape"], skew = coef(fitrSys)["skew"])

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

##################################### VaR/CoVaR evaluations ########################################
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

#Bivariate CoVaR: BTC-LTC
t_CoVaR_LTC_exc = evalBivariate(retX = rLTC, retY = rBTC, VaRX = VaR_rLTC,
                                CoVaR = t_CoVaR_1)
#Bivariate CoVaR: BTC-XMR
t_CoVaR_XMR_exc = evalBivariate(retX = rXMR, retY = rBTC, VaRX = VaR_rXMR,
                                CoVaR = t_CoVaR_2)
#Bivariate CoVaR: BTC-XRP
t_CoVaR_XRP_exc = evalBivariate(retX = rXRP, retY = rBTC, VaRX = VaR_rXRP,
                                CoVaR = t_CoVaR_3)
#SCoVaR
t_CoVaR_Sys_exc = evalBivariate(retX = rSys, retY = rBTC, VaRX = VaR_rSys,
                                CoVaR = t_CoVaR_4)
#MCoVaR
t_CoVaR_Multi_exc = evalMulti3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC, 
                                VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                MCoVaR = t_MCoVaR)
#VCoVaR
t_CoVaR_V_exc = evalVulnerability3X(retX1 = rLTC, retX2 = rXMR, retX3 = rXRP, retY = rBTC, 
                                    VaRX1 = VaR_rLTC, VaRX2 = VaR_rXMR, VaRX3 = VaR_rXRP, 
                                    VCoVaR = t_VCoVaR)

############################################# Plots ###############################################
VecToTS = function(vec, start = "2014-09-02", end = "2020-04-19"){
  inds <- seq(as.Date(start), as.Date(end), by = "day")
  return(zoo(as.vector(vec), inds))
}

time = decimal_date(time(VecToTS(rBTC)))
timeAT = c(which(time == 2015), which(time == 2016), which(time == 2017),
           which(time == 2018), which(time == 2019), which(time == 2020))
timeLAB = 2015:2020

############################################ BTC-LTC ##############################################
#create plot
plot.ts(cbind(rBTC, VaR_rBTC, t_CoVaR_1), plot.type = "single", col = c("blue3", "red3", "green3"),
        main = paste("Bivariate CoVaR of BTC|LTC - Realized violation rate: ", t_CoVaR_LTC_exc[[1]]), 
        ylim = c(-1.2, 0.25), ylab = "log-returns", xaxt = "n")
axis(1, at = timeAT, labels = timeLAB)

#cases, when BTC is below its VaR
points(VaR_BTC_exc, rep(-0.7, length(VaR_BTC_exc)), pch = 4)

#cases, when condition is fulfilled
points(VaR_LTC_exc, rep(-0.9, length(VaR_LTC_exc)), pch = 4)

#exceedences given the fulfilled condition
points(t_CoVaR_LTC_exc[[3]], rep(-1.1, t_CoVaR_LTC_exc[[2]]), pch = 4)

#further information
legend("bottom", horiz = TRUE, legend = c("log-return of BTC", "5% VaR of BTC", "5% CoVaR"), 
       fill = c("blue3", "red3", "green3"), bty = "n")
text(200, -0.65, expression("BTC"<="VaR(BTC):"), pos = 4, cex = 1)
text(200, -0.85, expression("LTC"<="VaR(LTC):"), pos = 4, cex = 1)
text(200, -1.05, expression(BTC<=CoVaR ~"|"~LTC <= paste(VaR(LTC), ":")), pos = 4, cex = 1)

############################################ BTC-XMR ##############################################
#create plot
plot.ts(cbind(rBTC, VaR_rBTC, t_CoVaR_2), plot.type = "single", col = c("blue3", "red3", "green3"),
        main = paste("Bivariate CoVaR of BTC|XMR - Realized violation rate: ", t_CoVaR_XMR_exc[[1]]), 
        ylim = c(-1.2, 0.25), ylab = "log-returns", xaxt = "n")
axis(1, at = timeAT, labels = timeLAB)

#cases, when BTC is below its VaR
points(VaR_BTC_exc, rep(-0.7, length(VaR_BTC_exc)), pch = 4)

#cases, when condition is fulfilled
points(VaR_XMR_exc, rep(-0.9, length(VaR_XMR_exc)), pch = 4)

#exceedences given the fulfilled condition
points(t_CoVaR_XMR_exc[[3]], rep(-1.1, t_CoVaR_XMR_exc[[2]]), pch = 4)

#further information
legend("bottom", horiz = TRUE, legend = c("log-return of BTC", "5% VaR of BTC", "5% CoVaR"), 
       fill = c("blue3", "red3", "green3"), bty = "n")
text(200, -0.65, expression("BTC"<="VaR(BTC):"), pos = 4, cex = 1)
text(200, -0.85, expression("XMR"<="VaR(XMR):"), pos = 4, cex = 1)
text(200, -1.05, expression(BTC<=CoVaR ~"|"~XMR <= paste(VaR(XMR), ":")), pos = 4, cex = 1)


############################################ BTC-XRP ##############################################
#create plot
plot.ts(cbind(rBTC, VaR_rBTC, t_CoVaR_3), plot.type = "single", col = c("blue3", "red3", "green3"),
        main = paste("Bivariate CoVaR of BTC|XRP - Realized violation rate: ", t_CoVaR_XRP_exc[[1]]), 
        ylim = c(-1.2, 0.25), ylab = "log-returns", xaxt = "n")
axis(1, at = timeAT, labels = timeLAB)

#cases, when BTC is below its VaR
points(VaR_BTC_exc, rep(-0.7, length(VaR_BTC_exc)), pch = 4)

#cases, when condition is fulfilled
points(VaR_XRP_exc, rep(-0.9, length(VaR_XRP_exc)), pch = 4)

#exceedences given the fulfilled condition
points(t_CoVaR_XRP_exc[[3]], rep(-1.1, t_CoVaR_XRP_exc[[2]]), pch = 4)

#further information
legend("bottom", horiz = TRUE, legend = c("log-return of BTC", "5% VaR of BTC", "5% CoVaR"), 
       fill = c("blue3", "red3", "green3"), bty = "n")
text(200, -0.65, expression("BTC"<="VaR(BTC):"), pos = 4, cex = 1)
text(200, -0.85, expression("XRP"<="VaR(XRP):"), pos = 4, cex = 1)
text(200, -1.05, expression(BTC<=CoVaR ~"|"~XRP <= paste(VaR(XRP), ":")), pos = 4, cex = 1)


########################################## BTC-System #############################################
#create plot
plot.ts(cbind(rBTC, VaR_rBTC, t_CoVaR_4), plot.type = "single", col = c("blue3", "red3", "green3"),
        main = paste("SCoVaR of BTC - Realized violation rate: ", t_CoVaR_Sys_exc[[1]]), 
        ylim = c(-1.2, 0.25), ylab = "log-returns", xaxt = "n")
axis(1, at = timeAT, labels = timeLAB)

#cases, when BTC is below its VaR
points(VaR_BTC_exc, rep(-0.7, length(VaR_BTC_exc)), pch = 4)

#cases, when condition is fulfilled
points(VaR_Sys_exc, rep(-0.9, length(VaR_Sys_exc)), pch = 4)

#exceedences given the fulfilled condition
points(t_CoVaR_Sys_exc[[3]], rep(-1.1, t_CoVaR_Sys_exc[[2]]), pch = 4)

#further information
legend("bottom", horiz = TRUE, legend = c("log-return of BTC", "5% VaR of BTC", "5% SCoVaR"), 
       fill = c("blue3", "red3", "green3"), bty = "n")
text(200, -0.65, expression("BTC"<="VaR(BTC):"), pos = 4, cex = 1)
text(200, -0.85, expression(sum(X[i])<="VaR("~sum(X[i])~"):"), pos = 4, cex = 1)
text(200, -1.05, expression(BTC<=SCoVaR ~"|"~sum(X[i])<="VaR("~sum(X[i])~"):"), pos = 4, cex = 1)


########################################### BTC-Multi ##############################################
Multi_exc = sort(intersect(VaR_LTC_exc, intersect(VaR_XMR_exc, VaR_XRP_exc)))

#create plot
plot.ts(cbind(rBTC, VaR_rBTC, t_MCoVaR), plot.type = "single", col = c("blue3", "red3", "green3"),
        main = paste("MCoVaR of BTC - Realized violation rate: ", t_CoVaR_Multi_exc[[1]]), 
        ylim = c(-1.75, 0.25), ylab = "log-returns", xaxt = "n")
axis(1, at = timeAT, labels = timeLAB)

#cases, when BTC is below its VaR
points(VaR_BTC_exc, rep(-1.2, length(VaR_BTC_exc)), pch = 4)

#cases, when condition is fulfilled
points(Multi_exc, rep(-1.4, length(Multi_exc)), pch = 4)

#exceedences given the fulfilled condition
points(t_CoVaR_Multi_exc[[3]], rep(-1.6, t_CoVaR_Multi_exc[[2]]), pch = 4)

#further information
legend("bottom", horiz = TRUE, legend = c("log-return of BTC", "5% VaR of BTC", "5% MCoVaR"), 
       fill = c("blue3", "red3", "green3"), bty = "n")
text(200, -1.125, expression("BTC"<="VaR(BTC):"), pos = 4, cex = 1)
text(200, -1.325, expression(symbol("\042")~i:~ X[i]<="VaR("~X[i]~"):"), pos = 4, cex = 1)
text(200, -1.525, expression(BTC<=MCoVaR ~"|"~symbol("\042")~i:~ X[i]<="VaR("~X[i]~"):"), pos = 4, cex = 1)


####################################### BTC-Vulnerability ##############################################
V_exc = sort(unique(c(VaR_LTC_exc, VaR_XMR_exc, VaR_XRP_exc)))

#create plot
plot.ts(cbind(rBTC, VaR_rBTC, t_VCoVaR), plot.type = "single", col = c("blue3", "red3", "green3"),
        main = paste("VCoVaR of BTC - Realized violation rate: ", t_CoVaR_V_exc[[1]]), 
        ylim = c(-1.15, 0.25), ylab = "log-returns", xaxt = "n")
axis(1, at = timeAT, labels = timeLAB)

#cases, when BTC is below its VaR
points(VaR_BTC_exc, rep(-0.6, length(VaR_BTC_exc)), pch = 4)

#cases, when condition is fulfilled
points(V_exc, rep(-0.8, length(V_exc)), pch = 4)

#exceedences given the fulfilled condition
points(t_CoVaR_V_exc[[3]], rep(-1.0, t_CoVaR_V_exc[[2]]), pch = 4)

#further information
legend("bottom", horiz = TRUE, legend = c("log-return of BTC", "5% VaR of BTC", "5% VCoVaR"), 
       fill = c("blue3", "red3", "green3"), bty = "n")
text(200, -0.55, expression("BTC"<="VaR(BTC):"), pos = 4, cex = 1)
text(200, -0.75, expression(symbol("\044")~i:~ X[i]<="VaR("~X[i]~"):"), pos = 4, cex = 1)
text(200, -0.95, expression(BTC<=VCoVaR ~"|"~symbol("\044")~i:~ X[i]<="VaR("~X[i]~"):"), pos = 4, cex = 1)

########################################## PANEL PLOT ##################################################

panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=4)
  txt <- paste0("R = ", r)
  text(0.5, 0.5, txt, cex = 1.5)
}

upper.panel<-function(x, y){
  points(x,y, pch = 19, col = "blue3")
}

# Create the plots
pairs(cbind("VaR(BTC)" = VaR_rBTC, "CoVaR|LTC" = t_CoVaR_1, "CoVaR|XMR" = t_CoVaR_2, 
            "CoVaR|XRP" = t_CoVaR_3, "SCoVaR" = t_CoVaR_4, "MCoVaR" = t_MCoVaR, "VCoVaR" = t_VCoVaR), 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
