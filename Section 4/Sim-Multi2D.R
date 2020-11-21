#Note: The calculations in this script take only a couple of seconds/minutes on an Intel(R) 
#      Xeon(R) Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("copula")

####################################### MULTI COVAR #######################################
get_MultiCoVaR = function(copula, theta, alpha = 0.05, beta = 0.05){
  if(copula == "gumbel"){
    cop = gumbelCopula(param = theta, dim = 3)
  }
  if(copula == "clayton"){
    cop = claytonCopula(param = theta, dim = 3)
  }
  
  MinF = function(v){
    return(pCopula(c(v, alpha, alpha), copula = cop) - pCopula(c(1, alpha, alpha), copula = cop)*beta)
  }
  tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                maxiter = 10000000)$root
  MCoVaR = qnorm(tmp)
  return(MCoVaR)
}

#Multi params
set.seed(1)
tau = seq(0.001, 0.97, length.out = 1000)

theta_Gu = iTau(copula = gumbelCopula(), tau = tau)
theta_Cl = iTau(copula = claytonCopula(), tau = tau)

### alpha = beta = 0.05 ###
#gumbel copula
M_gumbel = c()
for(i in seq_along(theta_Gu)){
  M_gumbel[i] = get_MultiCoVaR(copula = "gumbel", theta = theta_Gu[i], 
                               alpha = 0.05, beta = 0.05)
}

#clayton copula
M_clayton = c()
for(i in seq_along(theta_Cl)){
  M_clayton[i] = get_MultiCoVaR(copula = "clayton", theta = theta_Cl[i], 
                                alpha = 0.05, beta = 0.05)
}

#Comparison
plot(tau, M_clayton, type = "l", xlim = c(0, 1), ylab = "MCoVaR", col = "red3", lwd = 2,
     main = "Multi-CoVaR with standard normal margins", xlab = expression(tau))
lines(tau, M_gumbel, col = "blue3", lwd = 2)
legend("right", col = c("red3", "blue3"), bty = "n", density = c(0,0), lwd = 2,
       legend = c("clayton", "gumbel"), border = c(NA,NA), lty = c(1,1))
mtext(text = expression(paste(alpha, " = ", beta, " = 0.05")))
abline(h = qnorm(0.05*0.05), lwd = 2, col = "darkgrey", lty = 2)
abline(h = qnorm(0.05), lwd = 2, col = "darkgrey", lty = 2)

### alpha = beta = 0.01 ###
#gumbel copula
M_gumbel = c()
for(i in seq_along(theta_Gu)){
  M_gumbel[i] = get_MultiCoVaR(copula = "gumbel", theta = theta_Gu[i], 
                               alpha = 0.01, beta = 0.01)
}

#clayton copula
M_clayton = c()
for(i in seq_along(theta_Cl)){
  M_clayton[i] = get_MultiCoVaR(copula = "clayton", theta = theta_Cl[i], 
                                alpha = 0.01, beta = 0.01)
}

#Comparison
plot(tau, M_clayton, type = "l", xlim = c(0, 1), ylab = "MCoVaR", col = "red3", lwd = 2,
     main = "Multi-CoVaR with standard normal margins", xlab = expression(tau))
lines(tau, M_gumbel, col = "blue3", lwd = 2)
legend("right", col = c("red3", "blue3"), bty = "n", density = c(0,0), lwd = 2,
       legend = c("clayton", "gumbel"), border = c(NA,NA), lty = c(1,1))
mtext(text = expression(paste(alpha, " = ", beta, " = 0.01")))
abline(h = qnorm(0.01*0.01), lwd = 2, col = "darkgrey", lty = 2)
abline(h = qnorm(0.01), lwd = 2, col = "darkgrey", lty = 2)
