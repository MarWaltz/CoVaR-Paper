#Note: The calculations in this script take only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("tseries")
library("timeDate")

TestPipe = function(vec){
  vec = as.vector(vec)
  
  #Jarque-Bera test: H_0: data is normal distributed
  JB = jarque.bera.test(vec)$p.value
  
  #Ljung-Box test: H_0: no serial correlation
  LB = Box.test(vec, lag = 8, type = c("Ljung-Box"))$p.value
  
  #ADF-test: H_0: ts is non-stationary
  ADF = adf.test(vec, k = 8)$statistic
  
  #PP-test: H_0: ts is non-stationary
  PP = PP.test(vec)$statistic
  
  #KPSS-test: H_0: ts is (level) stationary
  KPSS = kpss.test(vec)$statistic
  
  return(round(c("Min" = min(vec),
                 "Mean" = mean(vec),
                 "Median" = median(vec),
                 "Max" = max(vec),
                 "Sd" = sd(vec),
                 "Kurtosis" = kurtosis(vec),
                 "Skewness" = skewness(vec),
                 "JB" = JB,
                 "LB" = LB,
                 "ADF" = ADF,
                 "PP" = PP,
                 "KPSS" = KPSS),4))
}

#load BTC, LTC, XMR, XRP from 'https://github.com/MarWaltz/Thesis/tree/master/Section 5/Data'
rBTC = diff(log(BTC))
rLTC = diff(log(LTC))
rXMR = diff(log(XMR))
rXRP = diff(log(XRP))
rSys = rLTC + rXMR + rXRP

#compute statistics and tests
res = cbind(TestPipe(rBTC), TestPipe(rLTC), TestPipe(rXMR), TestPipe(rXRP), TestPipe(rSys))
colnames(res) = c("BTC", "LTC", "XMR", "XRP", "System")
