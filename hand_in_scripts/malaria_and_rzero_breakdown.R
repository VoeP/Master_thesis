####malaria
##defining the functions and ranges

##more time steps so an equlibrium can be reached
time_steps <- seq(0,200,0.1)

x1<-seq(0.4, 20, by = 0.64)

sequenc<-seq(0, 3, by=0.3)

library(deSolve)
p <- list(alpha = 0.3, gamma=0.5, gammav=0, beta = 0.8, betav = 0.8, L=25, 
          lambdav=100, mu=0.1, muv=0.3, death=0)
print((p$beta*p$betav)/((p$gamma+p$mu)*(p$gammav+p$muv)))


gamma<-function(x){
  return(x/8)
}

beta<-function(x){
  max=4
  half=2
  return(((max*x)/(half+x)))
}

####defining the model (with density dependence)

f <- function(t, N, p){
  S1<-N[1]
  I1<-N[2]
  Sv1<-N[3]
  Iv1<-N[4]
  N1<-N[1]+N[2]
  Nv1<-N[3]+N[4]
  dS<- p$L-p$beta*S1*Iv1-(mort(N1)+0.1)*S1
  dI<- p$beta*S1*Iv1-(p$gamma + mort(N1)+0.1)*I1
  dSv<-p$lambdav-p$betav*Sv1*I1-p$muv*Sv1
  dIv<-p$betav*Sv1*I1-(p$muv+p$gammav)*Iv1
  return(list(c(dS, dI, dSv, dIv)))}



####defining more functions (population size, R0)
value_range<-seq(0, 3, 0.01)
RO<-function(x){
  return((beta(x)*beta(x)*(p$L/p$mu)*(p$lambdav/p$muv))/((p$mu+gamma(x))*p$muv))
}


S<-function(x){
  N0 <- c(10, 1, 10, 1)
  p <- list(alpha=0.3, gamma=gamma(x), gammav=0, beta = beta(x), betav = beta(x),
            L=25, lambdav=100, mu=0.1, muv=0.3, death=0)
  out<-ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
  return(c(tail(out[,2], n=1), tail(out[,3], n=1), tail(out[,4], n=1), tail(out[,5], n=1)))
}



####fitness function

fitness<-function(x1, x2, Sh, Sv, total2){
  return(beta(x2)*beta(x2)*Sh*Sv-((gamma(x2)+(mort(total2)+0.03))*p$muv))
}



p <- list(alpha = 1, gamma=0.5, gammav=0, beta = 0.8, betav = 0.8, L=25, 
          lambdav=100, mu=0.03, muv=0.3, death=0)

predictions<-list()
trues<-list()

##finding the virulences (first the true density dependent ESS
##and then later predictions)
for (h in 1:length(sequenc)){
  k=sequenc[h]
  mort<-function(x){
    
    return(k*x)
  }
  f <- function(t, N, p){
    S1<-N[1]
    I1<-N[2]
    Sv1<-N[3]
    Iv1<-N[4]
    N1<-N[1]+N[2]
    Nv1<-N[3]+N[4]
    dS<- p$L-p$beta*S1*Iv1-(mort(N1)+0.03)*S1
    dI<- p$beta*S1*Iv1-(p$gamma + mort(N1)+0.03)*I1
    dSv<-p$lambdav-p$betav*Sv1*I1-p$muv*Sv1
    dIv<-p$betav*Sv1*I1-(p$muv+p$gammav)*Iv1
    return(list(c(dS, dI, dSv, dIv)))}
  S<-function(x){
    N0 <- c(30, 1, 30, 1)
    p <- list(alpha=0.3, gamma=gamma(x), gammav=0, beta = beta(x), betav = beta(x),
              L=25, lambdav=100, mu=mort(32)+0.03, muv=0.3, death=0)
    out<-ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
    return(c(tail(out[,2], n=1), tail(out[,3], n=1), tail(out[,4], n=1), tail(out[,5], n=1)))
  }
  roi<-matrix(, nrow = length(x), ncol = length(x))
  x1<-seq(0.4, 20, by = 0.64)
  w1<-matrix(, nrow = length(x1), ncol = length(x1))
  for (i in 1:length(x1)){
    S_es<-S(x1[i])
    Sh1<-S_es[1]
    Si1<-S_es[2]
    Sv1<-S_es[3]
    Svi1<-S_es[4]
    total1<-S_es[1]+S_es[2]
    
    for (j in 1:length(x1)){
      S_es2<-S(x1[j])
      Sh2<-S_es2[1]
      Si2<-S_es2[2]
      Sv2<-S_es2[3]
      Svi2<-S_es2[4]
      total2<-S_es2[1]+S_es2[2]
      w1[i,j]=fitness(x1[i], x1[j], Sh1, Sv1, total2)
      print(paste(i, j))
      #w[i, j]<-(epsilon*((theta+d)/beta(x[i])))/(1+exp(-k*(x[j]-1.5)))-1-d
      if (RO(j)<1){
        break
      }
      if (RO(i)<1){
        break
      }
      print(paste(i, j))
      roi[i, j]<-total
      #    plot(1:41, (c((beta(x)*(p$L*(mort(total)+0.1)))/(gamma(x)+p$alpha+(mort(total)+0.1)))))
    }
  }
  clean_ro<-list()
  p <- list(alpha = 0.3, gamma=0.5, gammav=0, beta = 0.8, betav = 0.8, L=25, 
            lambdav=100, mu=0.03, muv=0.3, death=0)
  for (i in 1:length(x)){
    clean_ro[i]<-((beta(x[i])*beta(x[i])*(p$L/(mort(32)+0.03))*(p$lambdav/p$muv))/(((mort(32)+0.03)+gamma(x[i]))*p$muv))
  }
  #  plot(x, clean_ro)
  contour(x1, x1, w1)
  points_con<-contourLines(x1, x1, w1, levels=0)
  if(length(points_con)<2){
    next
  }

  ##The orientation of the contour lines was always the same, so the if-statements from the previous
  ##model are not necessary
  true_val<-max(points_con[[2]]$x)
  prediction<-x[which.max(clean_ro)]
  
  ##true value using the ESS 
  true<-min(points_con[[2]]$y)
  
  trues[h]<-true
  print(i)
}
predictions<-list()
for (i in 1:length(sequenc)){
  mort<-function(x){
    KN<-sequenc[i]*x
    return(KN)
  }
  f <- function(t, N, p){
    S1<-N[1]
    I1<-N[2]
    Sv1<-N[3]
    Iv1<-N[4]
    N1<-N[1]+N[2]
    Nv1<-N[3]+N[4]
    dS<- p$L-p$beta*S1*Iv1-p$mu*S1
    dI<- p$beta*S1*Iv1-(p$gamma + p$mu)*I1
    dSv<-p$lambdav-p$betav*Sv1*I1-p$muv*Sv1
    dIv<-p$betav*Sv1*I1-(p$muv+p$gammav)*Iv1
    return(list(c(dS, dI, dSv, dIv)))}
  
  S<-function(x){
    N0 <- c(30, 0, 30, 0)
    p <- list(alpha=0.3, gamma=gamma(x), gammav=0, beta = beta(x), betav = beta(x),
              L=25, lambdav=100, mu=mort(30)+0.03, muv=0.3, death=0)
    out<-ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
    return(c(tail(out[,2], n=1), tail(out[,3], n=1), tail(out[,4], n=1), tail(out[,5], n=1)))
  }
  S2<-function(x){
    N0 <- c(30, 1, 30, 1)
    p <- list(alpha<-0.3, gamma=gamma(x), gammav=0, beta = beta(x), betav = beta(x),
              L=25, lambdav=100, mu=mort(31)+0.03, muv=0.3, death=0)
    out<-ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
    return(c(tail(out[,2], n=1), tail(out[,3], n=1), tail(out[,4], n=1), tail(out[,5], n=1)))
  }
  
  ros<-list()
  ##getting necessary population sizes for the R0:formula
  for (j in 1:length(x1)){
    print(j/length(x1))
    S_es<-S(0)
    print(S_es[1])
    total<-S_es[1]+S_es[2]
    S_es2<-S2(x1[j])
    print(S_es2)
    total<-S_es[1]+S_es[2]
    ##finding the R0:s using formula
    ros[j]<-(beta(x1[j])*beta(x1[j])*(S_es[1])*(S_es[3]))/((mort(31)+0.03+gamma(x1[j]))*p$muv)
  }
  ##prediction using the maximal R0
  predictions[i]<-x1[which.max(ros)]
}

plot(sequenc[1:length(predictions)], trues[1:length(predictions)], type="l", col="red", xlim = c(0, 0.5), ylim = c(2,10), xlab="relationship variable", ylab="virulence")
lines(sequenc[1:length(predictions)], predictions[1:length(predictions)], col="blue")



##moving on to the R0:s

f <- function(t, N, p){
  S1<-N[1]
  I1<-N[2]
  Sv1<-N[3]
  Iv1<-N[4]
  N1<-N[1]+N[2]
  Nv1<-N[3]+N[4]
  dS<- p$L-p$beta*S1*Iv1-(mort(N1)+0.1)*S1
  dI<- p$beta*S1*Iv1-(p$gamma + mort(N1)+0.1)*I1
  dSv<-p$lambdav-p$betav*Sv1*I1-p$muv*Sv1
  dIv<-p$betav*Sv1*I1-(p$muv+p$gammav)*Iv1
  return(list(c(dS, dI, dSv, dIv)))}



trues_ros<-list()


p <- list(alpha = 1, gamma=0.5, gammav=0, beta = 0.8, betav = 0.8, L=25, 
          lambdav=100, mu=0.03, muv=0.3, death=0)



S<-function(x){
  N0 <- c(30, 1, 30, 1)
  p <- list(alpha=0.3, gamma=gamma(x), gammav=0, beta = beta(x), betav = beta(x),
            L=25, lambdav=100, mu=0.1, muv=0.3, death=0)
  out<-ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
  return(c(tail(out[,2], n=1), tail(out[,3], n=1), tail(out[,4], n=1), tail(out[,5], n=1)))
}



for (i in 1:length(trues)){
  mort<-function(x){
    KN<-sequenc[i]*x
    return(KN)
  }
  S_es<-S(0)
  print(S_es[1])
  total<-S_es[1]+S_es[2]
  S_es2<-S(trues[[i]])
  print(S_es2)
  total<-S_es[1]+S_es[2]
  #    ros[j]<-(S_es[1]*S_es[3])/(S_es2[1]*S_es2[3])
  ##so: question: is it the mort(S_es2 population), or is it mort(31)? Makes a pretty big difference actually, but if it is the R0, then it shouldbethe mort(31)Ithink
  trues_ros[i]<-(beta(trues[[i]])*beta(trues[[i]])*(S_es[1])*(S_es[3]))/((mort(31)+0.03+gamma(trues[[i]]))*p$muv)
  
}



predictions_ros<-list()
for (i in 1:length(predictions)){
  mort<-function(x){
    KN<-sequenc[i]*x
    return(KN)
  }
  S_es<-S(0)
  print(S_es[1])
  total<-S_es[1]+S_es[2]
  S_es2<-S(predictions[[i]])
  print(S_es2)
  total<-S_es[1]+S_es[2]
  #    ros[j]<-(S_es[1]*S_es[3])/(S_es2[1]*S_es2[3])
  predictions_ros[i]<-(beta(predictions[[i]])*beta(predictions[[i]])*(S_es[1])*(S_es[3]))/((mort(31)+0.03+gamma(predictions[[i]]))*p$muv)
  
}




plot(sequenc, predictions_ros, xlim=c(0, 0.5),  type="l", col="blue")
lines(sequenc, trues_ros, col="red")

