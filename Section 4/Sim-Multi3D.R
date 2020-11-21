#Note: The calculations in this file take approximately 0.5h on an Intel(R) Xeon(R) 
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("HAC")
library("plotly")

############################# 3D Plot HAC MULTI COVAR ##########################
get_MultiCoVaR3D = function(copula, theta1, theta2, alpha = 0.05, beta = 0.05){
  if(theta1 > theta2){
    return(NA)}
  if(copula == "gumbel"){
    hac = hac(type = 1, tree = list("v", list("u1", "u2", theta2), theta1))
  }
  if(copula == "clayton"){
    hac = hac(type = 3, tree = list("v", list("u1", "u2", theta2), theta1))
  }
  
  MinF = function(v){
    return(pHAC(X = c("v" = v, "u1" = alpha, "u2" = alpha), hac = hac)- 
             pHAC(X = c("v" = 1, "u1" = alpha, "u2" = alpha), hac = hac)*beta)
  }
  tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
  
  MCoVaR = qnorm(tmp)
  return(MCoVaR)
}

length = 600; set.seed(1)
tau = seq(0.001,0.98,length.out = length)

#Gumbel
theta_Gu = iTau(copula = gumbelCopula(), tau = tau)
z_Gu = matrix(NA, nrow = length, ncol = length)

for(i in 1:nrow(z_Gu)){
  for(j in 1:ncol(z_Gu)){
    z_Gu[i,j] = get_MultiCoVaR3D(copula = "gumbel", theta1 = theta_Gu[i],
                                 theta2 = theta_Gu[j], alpha = 0.05, beta = 0.05)
  }
}

#Clayton
theta_Cl = iTau(copula = claytonCopula(), tau = tau)
z_Cl = matrix(NA, nrow = length, ncol = length)

for(i in 1:nrow(z_Cl)){
  for(j in 1:ncol(z_Cl)){
    z_Cl[i,j] = get_MultiCoVaR3D(copula = "clayton", theta1 = theta_Cl[i],
                                 theta2 = theta_Cl[j], alpha = 0.05, beta = 0.05)
  }
}

#Surface plot
plot_ly()%>%
  add_surface(x=~tau, y=~tau, z =~t(z_Gu),
              opacity = 0.95,
              showscale = T,
              colorscale = "Blues",
              colorbar = list(title = "gumbel"),
              contours = list(z = list(show = F, 
                                       usecolormap = TRUE,
                                       project = list(z=TRUE))
              ))%>%
  add_surface(x=~tau, y=~tau, z =~t(z_Cl),
              opacity = 1,
              showscale = T,
              colorscale = "Reds",
              colorbar = list(title = "clayton"),
              reversescale  = T,
              contours = list(z = list(show = F, 
                                       usecolormap = TRUE,
                                       project = list(z=TRUE))
              ))%>% 
  layout(scene = list(camera = list(eye = list(x = 2.25, y = -1.75, z = 0.75)),
                      xaxis = list(title = "tau1"),
                      yaxis = list(title = "tau2"),
                      zaxis = list(title = "MCoVaR")))
