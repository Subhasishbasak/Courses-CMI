## What is the average blood pressure reading of a sub-population of 20 adult men?


ybar<-128
s<-7.67
n<-20
sim.size<-10000
set.seed(2928)

## Analysis with Flat prior or Non-informative Improper Prior
### Simulate mu
t<-rt(sim.size,df=n-1)
mu<-ybar+t*s/sqrt(n)
cat("Summary of mu :","\n")
summary(mu)
cat("Variance of mu :","\n")
var(mu)
cat("Sd of mu :","\n")
sd(mu)
cat("95% CI of mu :","\n")
quantile(mu,prob=c(0.025,0.975))

## Application: Expert's Opinion with Conjugate Prior
mu_0=120
sigma_0=5
k_0=1
nu_0=1

mu_n<-(k_0/(k_0+n))*mu_0+(n/(k_0+n))*ybar
k_n<-k_0+n
nu_n<-nu_0+n
sigma_n.sq<-(nu_0*sigma_0^2+(n-1)*s^2+(k_0*n/(k_0+n))*(ybar-mu_0)^2)/nu_n
sigma_n<-sqrt(sigma_n.sq)
sim.size<-10000
set.seed(2928)
### Simulate mu
t<-rt(sim.size,df=nu_n)
mu<-ybar+t*sigma_n/sqrt(k_n)
cat("Summary of mu :","\n")
summary(mu)
cat("Variance of mu :","\n")
var(mu)
cat("Sd of mu :","\n")
sd(mu)
cat("95% CI of mu :","\n")
quantile(mu,prob=c(0.025,0.975))

## Now try n = 4,6,8,10,12,16,20
## For each choices of n note down the 95% CI for both flat prior and conjugate prior method.
