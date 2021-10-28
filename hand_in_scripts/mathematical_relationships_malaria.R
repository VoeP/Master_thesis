####finding mathematical relationships
time_steps <- seq(0,150,0.5)

library(deSolve)
p <- list(alpha = 1, gamma=0.5, gammav=0, beta = 0.8, betav = 0.8, L=25, 
          lambdav=25, mu=0.01, muv=0.01, death=0)
print((p$beta*p$betav)/((p$gamma+p$mu)*(p$gammav+p$muv)))
value_range<-range(0, 3, 0.1)

gamma<-function(x){
  return( x)
}

beta<-function(x){
  epsilon=10
  k=8
  x0 = 0.5
  return((epsilon/(1+exp(-k*(x-x0))))/2)
}
for (i in value_range){
  print(beta(i))
}

beta_2<-function(x){
  epsilon=10
  k=2
  x0 = 0.5
  return((1))
}
for (i in value_range){
  print(beta(i))
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

value_range<-seq(0, 3, 0.1)
RO<-function(x){
  return((beta(x)*beta(x)*(p$L/p$mu)*(p$lambdav/p$muv))/((p$mu+gamma(x))*p$muv))
}


S<-function(x){
  N0 <- c(20, 5, 20, 5)
  p <- list(alpha = 1, gamma=gamma(x), gammav=0, beta = beta(x), betav = beta(x),
            L=25, lambdav=25, mu=0.01, muv=0.01, death=0)
  out<-ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
  return(c(tail(out[,2], n=1), tail(out[,4], n=1)))
}


I<-function(x){
  N0 <- c(20, 5, 20, 5)
  p <- list(alpha = 1, gamma=gamma(x), gammav=0, beta = beta(x), betav = beta(x),
            L=500, lambdav=500, mu=0.1, muv=0.1, death=0)
  out<-ode(y = N0, times = time_steps, func = f, parms = p, method = c("ode45"))
  return(c(tail(out[,3], n=1), tail(out[,5], n=1)))
}


Sv_s<-list()
Sh_s<-list()
for (i in 1:length(value_range)){
  S_es<-S(value_range[i])
  Sv_s[i]<-S_es[2]
  Sh_s[i]<-S_es[1]
}


ROs<-list()
for (i in 1:length(value_range)){
  ROs[i]<- RO(value_range[i])
}

plot(value_range, Sv_s, type="l")
plot(value_range, Sh_s, type="l")

par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(Sh_s, Sv_s, type="l", col="blue") # first plot
legend("topright",legend=c("vector population against", "susceptible human population"),
       text.col=c("red", "red"))

ShSv<-list()
for (i in 1:length(Sv_s)){
  ShSv[i]<-Sv_s[[i]]*Sh_s[[i]]
}

par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(value_range, Sv_s, type="l", col="blue") # first plot

par(new = TRUE)
plot(value_range, Sh_s, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="red")
mtext("Susceptible humans", side=4, line=3)

par(new = TRUE)
plot(value_range, ShSv, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="purple")



##here, I was just looking at the R0 and the proxy

ROx<-function(x){
  return(((beta(x)*beta(x)*(p$L/p$mu)*(p$lambdav/p$muv))/((p$mu+gamma(x))*p$muv)-1))*((gamma(x)+p$mu)*p$muv)
}
proxyvalues<-list()
for (i in 1:length(value_range)){
  proxyvalues[i]<-ROx(value_range[i])
}


##This plot shows that with the peak of R0, the product of the populations is minimised.
##Similar to what is shown in figure 3.
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(value_range, ROs, type="l", col="blue") # first plot
par(new = TRUE)
plot(value_range, ShSv, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="purple")
abline(v=0.8)




rO_division<-list()
S_division<-list()
for (i in 1:length(proxyvalues)){
  rO_division[i]<-proxyvalues[[i]]/proxyvalues[[i+1]]
}

for (i in 1:length(ShSv)){
  S_division[i]<-ShSv[[i+1]]/ShSv[[i]]
}

##This plot shows the mathematical relationship shown in figure 4 of the thesis. The first index is skipped,
##since there are some indexing issues
plot(value_range[2:30], rO_division[2:30], xlab = "virulence", ylab = "R0 mutant/ R0 resident", type="l", col="blue")
plot(value_range[2:30], S_division[2:30],  xlab = "virulence", ylab = "Susceptible resident/ susceptible mutant", type="l", col="red")
