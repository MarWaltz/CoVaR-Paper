#Note: The calculations in this script take only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

########################################## PLOT DATA #############################################
library("zoo")
library("lubridate")

VecToTS = function(vec, start = "2014-09-02", end = "2020-04-19"){
  inds <- seq(as.Date(start), as.Date(end), by = "day")
  return(zoo(as.vector(vec), inds))
}

MaxMinScale = function(ts){
  return((ts - min(ts))/(max(ts)-min(ts)))
}

#load BTC, LTC, XMR, XRP from 'https://github.com/MarWaltz/Thesis/tree/master/Section 5/Data'
rBTC = diff(log(BTC))
rLTC = diff(log(LTC))
rXMR = diff(log(XMR))
rXRP = diff(log(XRP))

#plot
par(mfrow = c(1,2))
timeLAB = 2015:2020

time = decimal_date(time(VecToTS(BTC, start = "2014-09-01")))
timeAT = c(which(time == 2015), which(time == 2016), which(time == 2017),
           which(time == 2018), which(time == 2019), which(time == 2020))
plot.ts(cbind(MaxMinScale(BTC), MaxMinScale(LTC), MaxMinScale(XMR), MaxMinScale(XRP)), xaxt = "n",
        plot.type = "single", col = c("blue3", "green3", "red3", "orange"), ylab = "Standardized Prices")
axis(1, at = timeAT, labels = timeLAB)
legend("topleft", col = c("blue3", "green3", "red3", "orange3"), bty = "n",
       legend = c("BTC", "LTC", "XMR", "XRP"), border = rep(NA, 4), lty = rep(1, 4), 
       density = rep(0, 4), lwd = 2)

time = decimal_date(time(VecToTS(rBTC, start = "2014-09-02")))
timeAT = c(which(time == 2015), which(time == 2016), which(time == 2017),
           which(time == 2018), which(time == 2019), which(time == 2020))
plot.ts(cbind(rBTC, rLTC, rXMR, rXRP), plot.type = "single", col = c("blue3", "green3", "red3", "orange3"),
     ylab = "log-returns", xaxt = "n")
axis(1, at = timeAT, labels = timeLAB)
legend("topleft", col = c("blue3", "green3", "red3", "orange3"), bty = "n",
       legend = c("BTC", "LTC", "XMR", "XRP"), border = rep(NA, 4), lty = rep(1, 4), 
       density = rep(0, 4), lwd = 2)
