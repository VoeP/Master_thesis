####malaria with numbers instead of proproptions
library(deSolve)


Sh<-function(x){
  N0 <- c(20, 5, 20, 5)
  p <- list(alpha = 1, gamma=gamma(x), gammav=0, beta = beta(x), betav = beta(x),
            L=12.5, lambdav=25, mu=0.1, muv=0.3, death=0)
  out<-ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
  return(tail(out[,2], n=1))
}
Sv<-function(x){
  N0 <- c(20, 5, 20, 5)
  p <- list(alpha = 1, gamma=gamma(x), gammav=0, beta = beta(x), betav = beta(x),
            L=12.5, lambdav=25, mu=0.1, muv=0.3, death=0)
  out<-ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
  return(tail(out[,4], n=1))
}


##the malaria-model
f4 <- function(t, N, p){
  S1<-N[1]
  I1<-N[2]
  Sv1<-N[3]
  Iv1<-N[4]
  Ihm<-N[5]
  Ivm<-N[6]
  N1<-N[1]+N[2]+N[5]
  Nv1<-N[3]+N[4]+N[6]
  #Nhm<-N[1]+N[5]
  #Nvm<-N[3]+N[6]
  #dS<- p$L-beta(p$v1)*S1*(Iv1/Nv1)-p$mu*S1-beta(p$v2)*S1*(Ivm/Nv1)
  #dI<- beta(p$v1)*S1*(Iv1/Nv1)-(gamma(p$v1) + p$mu)*I1
  #dSv<-p$lambdav-beta(p$v1)*Sv1*(I1/N1)-p$muv*Sv1-beta(p$v2)*Sv1*(Ihm/N1)
  #dIv<-beta(p$v1)*Sv1*(I1/N1)-(p$muv+gamma(p$v1))*Iv1
  #dIhm<- beta(p$v2)*S1*(Ivm/Nv1)-(gamma(p$v2) + p$muHm)*Ihm
  #dIvm<-beta(p$v2)*Sv1*(Ihm/N1)-(p$muVm+gamma(p$v2))*Ivm
  dS<- p$L-p$beta*S1*Iv1-p$mu*S1-p$betaHm*S1*Ivm
  dI<- p$beta*S1*Iv1-(p$gamma + p$mu)*I1
  dSv<-p$lambdav-p$betav*Sv1*I1-p$muv*Sv1-p$betaVm*Sv1*Ihm
  dIv<-p$betav*Sv1*I1-(p$muv+p$gammav)*Iv1
  dIhm<- p$betaHm*S1*Ivm-(p$gammaHm + p$muHm)*Ihm
  dIvm<-p$betaVm*Sv1*Ihm-(p$muVm+p$gammaVm)*Ivm
  return(list(c(dS, dI, dSv, dIv, dIhm, dIvm)))}

gamma<-function(x){
  return( x)
}

beta<-function(x){
  epsilon=10
  k=5
  x0 = 0.5
  return((epsilon/(1+exp(-k*(x-x0))))/2)
}



library(deSolve)
value_range<-seq(0, 1.5, 0.1)
matches<-list()
R_zeros<-list()
current<-0
N0 <- c(5, 15, 5, 15, 1, 1)
time_steps <- seq(0,1500,0.05)
for (i in 1:16){
  new_vir<-value_range[i]
  newbeta<-beta(new_vir)
  old_beta<-beta(current)
  new_gamma<-gamma(new_vir)
  old_gamma<-gamma(current)
  print(i)
  q <- list(alpha = 1, gamma=old_gamma, gammav=0, beta = old_beta, betav = old_beta, L=12.5, lambdav=25, mu=0.1, muv=0.3,
            death=0, muVm=0.3, muHm=0.1, gammaVm=0, gammaHm=new_gamma, betaVm=newbeta, betaHm=newbeta, 
            v1=current, v2=new_vir)
  
  eq<-ode(y = N0, times = time_steps, func = f4, parms = q, method = c("ode45"))
  print("condition")
  #if (tail(eq[,5], n=1)<2 && tail(eq[,7], n=1)>2 && tail(eq[,2], n=1)>1 && tail(eq[,6], n=1)>1){
  if (tail(eq[,5], n=1)<2 && tail(eq[,7], n=1)>2 ){
    print("Match")
    current<-new_vir
    matches[i]<-1
    N0<-c(tail(eq[,2], n=1), tail(eq[,6], n=1), tail(eq[,4], n=1), tail(eq[,7], n=1), 1, 1)
  }  else{
    matches[i]<-0
  }
  R_zeros[i]<-R0_v2(value_range[i])
  plot(NA, xlim=c(0,400), ylim=c(0, 30), ylab = "")
  lines(x = eq[,1], y = eq[,2], col='blue')
  lines(x = eq[,1], y = eq[,3], col='red')
  lines(x = eq[,1], y = eq[,4], col='yellow')
  lines(x = eq[,1], y = eq[,5], col='green')
  lines(x = eq[,1], y = eq[,6], col=645)
  lines(x = eq[,1], y = eq[,7], col='purple')
  abline(h=1)
  abline(h=0)
  print("current value", value_range[i])
}

##here I show that the value at which the new mutant will be unable to replace the resident is
##also the value at which the R0 is maximized, and thus where the ESS would be
plot(value_range, matches, type = "l")
plot(value_range, R_zeros, type = "l")
abline(v=0.9502)