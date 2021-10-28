####SIR model with recovery included


##note: the system gets run twice basically. This could have been rewritten
##so that everything gets done at the same time, but I didn't want to risk confucing myself 
##and losing a lot of time

library(deSolve)

x1<-seq(0, 3, 0.1)
w1<-seq(0, 3, 0.1)
gamma<-function(x){
  return( x)
}

mort<-function(x){
  KN<-0.1*x
  return(KN)
}

#beta<-function(x){
#   max=4
#  half=2
#  return(((max*x)/(half+x)))
#}


beta<-function(x){
  epsilon=1
  k=5
  x0 = 0.5
  return((epsilon/(1+exp(-k*(x-x0))))/2)
}

#beta<-function(x){
#  power<-2
#  return(x^power)
#}


##plotting the beta-function

list_of_betas<-list()
for (i in 1:length(x)){
  list_of_betas[i]<-beta(x[i])
}
plot(x, list_of_betas)


time_steps <- seq(0,50,0.05)
p <- list(beta = 0.1, gamma = 0.05, mu = (0.1+mort(25)), L=12, alpha=0.1, delta=0.05)
#N0 <- c(20, 5, 0, 25)

alpha<-function(x){
  return(0.1)
}




####finding the virulences and R0:s

predictions<-list()
trues<-list()

sequenc<-seq(0, 0.5, by=0.0075)

for (h in 1:length(sequenc)){
  k=sequenc[h]
  mort<-function(x){
    
    return(k*x)
  }
  ##system with density dependence incorporated
  f <- function(t, N, p){
    S1<-N[1]
    I1<-N[2]
    R1<-N[3]
    Ntot <-S1 + I1 + R1
    dS<- p$L-p$beta*S1*I1-(0.1+k*Ntot)*S1+p$delta*R1
    dI<- p$beta*S1*I1-p$gamma*I1-(0.1+k*Ntot)*I1-p$alpha*I1
    dR<-p$alpha*I1-(0.1+k*Ntot)*R1-p$delta*R1
    return(list(c(dS, dI, dR)))}
  
  p <- list(beta = 0.1, gamma = 0.05, mu = (0.1+mort(25)), L=12, alpha=0.1, delta=0.05)
  S_s<-function(x){
    N0 <- c(20, 5, 1)
    p <- list(alpha = 0.1, beta = beta(x), 
              gamma = gamma(x), mu = (0.1+mort(25)), L=12, delta=0.05)
    eq <- ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
    return(c(tail(eq[,2], n = 1), tail(eq[,3], n = 1), tail(eq[,4], n = 1)))
  }
  x<-seq(0, 3, 0.05)
  w<-matrix(, nrow = length(x), ncol = length(x))
  roi<-matrix(, nrow = length(x), ncol = length(x))
  for (i in 1:length(x)){
    S_pop<-S_s(x[i])
    total<-S_pop[1]+S_pop[2]+S_pop[3]
    for (j in 1:length(x)){
      S_pop2<-S_s(x[j])
      total2<-S_pop2[1]+S_pop2[2]+S_pop2[3]
      w[i,j]=(beta(x[j])*S_pop[1])-(gamma(x[j])+p$alpha+(mort(total2)+0.1))
      #w[i, j]<-(epsilon*((theta+d)/beta(x[i])))/(1+exp(-k*(x[j]-1.5)))-1-d
      print(paste(i, j))
      roi[i, j]<-total
      #    plot(1:41, (c((beta(x)*(p$L*(mort(total)+0.1)))/(gamma(x)+p$alpha+(mort(total)+0.1)))))
    }
  }
  p <- list(beta = 0.1, gamma = 0.05, mu = (0.1+mort(25)), L=12, alpha=0.1, delta=0.05)
  clean_ro<-list()
  for (i in 1:length(x)){
    S_pop<-S_s(x[i])
    total<-S_pop[1]+S_pop[2]+S_pop[3]
    clean_ro[i]<-((beta(x[i])*((0.1+sqrt(1+4*k*p$L))/2)/(gamma(x[i])+p$alpha+p$mu)))
  }
  plot(x, clean_ro)
  
  ##finding the ESS using the contour lines and their orientation.
  points_con<-contourLines(x, x, w, levels=0)
  contour(x, x, w, levels = 0, drawlabels = TRUE)
  lines(x=points_con[[1]]$x, y=points_con[[1]]$y, col="red")
  lines(x=points_con[[2]]$x, y=points_con[[2]]$y, col="green")
  plot(points_con[[1]]$x, points_con[[1]]$y, type="l")
  plot(points_con[[2]]$x, points_con[[2]]$y, type="l")
  true_val<-max(points_con[[2]]$x)
  prediction<-x[which.max(clean_ro)]
  predictions[h]<-prediction
  ##the if-statement will find the ESS based on how R identifies the contour-lines
  if (max(points_con[[1]]$y>1.8)){
    true<-min(points_con[[2]]$x)
  } else{
    if (max(points_con[[1]]$x>1.8)){
      true<-max(points_con[[1]]$y)
    }
    
  }
  trues[h]<-true
  prediction<-x[which.max(clean_ro)]
  predictions[h]<-prediction
}

##plotting the virulences
plot(sequenc[1:length(trues)], predictions[1:length(trues)], type="l", col="blue")
plot(sequenc[1:length(trues)], trues[1:length(trues)], type="l", col="blue")
plot(predictions[1:length(trues)], trues[1:length(trues)], type="l", col="blue")


plot(sequenc[3:length(trues)], predictions[3:length(trues)], ylim=c(0.8, 1.15), type="l", col="blue")
lines(sequenc[3:length(trues)], trues[3:length(trues)], col="red")



abline(v=0.55)
abline(h=0.55)

differences<-list()
for (i in 1:length(trues)){
  differences[i]<-predictions[[i]]-trues[[i]]
}
plot(sequenc[1:length(trues)], differences[1:length(trues)], type="l")





##getting the R0:s
sequenc<-seq(0, 0.3, by=0.02)
predictions2<-list()
for (h in 1:length(sequenc)){
  k=sequenc[h]
  mort<-function(x){
    
    return(k*x)
  }
  ##again, the system with incorporated density-dependence
  f <- function(t, N, p){
    S1<-N[1]
    I1<-N[2]
    R1<-N[3]
    Ntot <-S1 + I1 + R1
    dS<- p$L-p$beta*S1*I1-(0.1+k*Ntot)*S1+p$delta*R1
    dI<- p$beta*S1*I1-p$gamma*I1-(0.1+k*Ntot)*I1-p$alpha*I1
    dR<-p$alpha*I1-(0.1+k*Ntot)*R1-p$delta*R1
    return(list(c(dS, dI, dR)))}
  
  p <- list(beta = 0.1, gamma = 0.05, mu = (0.1+mort(25)), L=12, alpha=0.1, delta=0.05)
  S_s<-function(x){
    N0 <- c(20, 5, 1)
    p <- list(alpha = 0.1, beta = beta(x), 
              gamma = gamma(x), mu = (0.1+mort(25)), L=12, delta=0.05)
    eq <- ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
    return(c(tail(eq[,2], n = 1), tail(eq[,3], n = 1), tail(eq[,4], n = 1)))
  }
  x<-seq(0, 3, 0.05)
  w<-matrix(, nrow = length(x), ncol = length(x))
  roi<-matrix(, nrow = length(x), ncol = length(x))
  p <- list(beta = 0.1, gamma = 0.05, mu = (0.1+mort(25)), L=12, alpha=0.1, delta=0.05)
  clean_ro<-list()
  for (i in 1:length(x)){
    S_de<-S_s(0)
    Ses<-S_s(x[i])
    ##finding the R0 based on Shat/S*
    clean_ro[i]<-(S_de[1]/Ses[1])
  }
  prediction<-x[which.max(clean_ro)]
  predictions2[h]<-prediction
  
}
plot(1:16, predictions2)

ro_list_trues<-list()
ro_list_predictions<-list()
for (i in 1:30){
  k=sequenc[i]
  mort<-function(x){
    
    return(k*x)
  }
  ##the system with incorporated density-dependence
  f <- function(t, N, p){
    S1<-N[1]
    I1<-N[2]
    R1<-N[3]
    Ntot <-S1 + I1 + R1
    dS<- p$L-p$beta*S1*I1-(0.1+k*Ntot)*S1+p$delta*R1
    dI<- p$beta*S1*I1-p$gamma*I1-(0.1+k*Ntot)*I1-p$alpha*I1
    dR<-p$alpha*I1-(0.1+k*Ntot)*R1-p$delta*R1
    return(list(c(dS, dI, dR)))}
  
  p <- list(beta = 0.1, gamma = 0.05, mu = (0.1+mort(25)), L=12, alpha=0.1, delta=0.05)
  S_s<-function(x){
    N0 <- c(20, 5, 1)
    p <- list(alpha = 0.1, beta = beta(x), 
              gamma = gamma(x), mu = (0.1+mort(25)), L=12, delta=0.05)
    eq <- ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
    return(c(tail(eq[,2], n = 1), tail(eq[,3], n = 1), tail(eq[,4], n = 1)))
  }
  S_de<-S_s(0)
  #  Ses<-S_s(trues[[i]])
  ro_list_trues[i]<-((beta(trues[[i]])*((S_de[1]))/(gamma(trues[[i]])+p$alpha+p$mu)))
}
for (i in 1:30){
  k=sequenc[i]
  mort<-function(x){
    
    return(k*x)
  }
  f <- function(t, N, p){
    S1<-N[1]
    I1<-N[2]
    R1<-N[3]
    Ntot <-S1 + I1 + R1
    dS<- p$L-p$beta*S1*I1-(0.1+k*Ntot)*S1+p$delta*R1
    dI<- p$beta*S1*I1-p$gamma*I1-(0.1+k*Ntot)*I1-p$alpha*I1
    dR<-p$alpha*I1-(0.1+k*Ntot)*R1-p$delta*R1
    return(list(c(dS, dI, dR)))}
  
  p <- list(beta = 0.1, gamma = 0.05, mu = (0.1+mort(25)), L=12, alpha=0.1, delta=0.05)
  S_s<-function(x){
    N0 <- c(20, 5, 1)
    p <- list(alpha = 0.1, beta = beta(x), 
              gamma = gamma(x), mu = (0.1+mort(25)), L=12, delta=0.05)
    eq <- ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
    return(c(tail(eq[,2], n = 1), tail(eq[,3], n = 1), tail(eq[,4], n = 1)))
  }
  S_de<-S_s(0)
  Ses<-S_s(predictions[[i]])
  ro_list_predictions[i]<-(S_de[1]/Ses[1])
}

plot(sequenc[1:length(ro_list_trues)], ro_list_trues, col="blue", type="l")
lines(sequenc[1:length(ro_list_trues)], ro_list_predictions, col="red")

plot(sequenc[1:16], ro_list_trues[1:16], col="blue", type="l")
lines(sequenc[1:16], ro_list_predictions[1:16], col="red")

plot(sequenc[1:16], ro_list_trues[1:16], col="red", type="l", xlim=c(0, 0.3), ylim = c(0,5.5))
lines(sequenc[1:16], ro_list_predictions[1:16], col="blue")

